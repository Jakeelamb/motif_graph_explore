/// K-mer counting with parallel processing, working directly on mmap bytes.
///
/// Counts k-mer occurrences across genome sequences without copying sequence
/// data — newlines are skipped inline during the sliding window scan.

use crate::encode::{self, KmerU64};
use crate::fasta::FastaReader;
use rayon::prelude::*;
use rustc_hash::FxHashMap;

/// The vocabulary: k-mers at various lengths that passed filtering.
pub struct Vocabulary {
    entries: FxHashMap<(u64, u8), VocabEntry>,
}

#[derive(Clone, Debug)]
pub struct VocabEntry {
    pub kmer: KmerU64,
    pub count: u32,
}

/// Base frequencies computed from k=1 counts.
/// Used for composition-aware significance testing.
pub struct BaseFreqs {
    /// Frequency of each base, indexed by 2-bit encoding: A=0, C=1, G=2, T=3
    pub freqs: [f64; 4],
}

impl BaseFreqs {
    /// Compute base frequencies from k=1 raw counts.
    /// Encoding: A=0, C=1, G=2, T=3.
    pub fn from_counts(counts_1: &FxHashMap<u64, u32>) -> Self {
        let mut raw = [0u64; 4];
        for (&encoded, &count) in counts_1 {
            if encoded < 4 {
                raw[encoded as usize] += count as u64;
            }
        }
        let total = raw.iter().sum::<u64>() as f64;
        let freqs = if total > 0.0 {
            [
                raw[0] as f64 / total,
                raw[1] as f64 / total,
                raw[2] as f64 / total,
                raw[3] as f64 / total,
            ]
        } else {
            [0.25; 4]
        };
        BaseFreqs { freqs }
    }

    /// Expected count under independence model: λ = (L - k + 1) × ∏ p(base_i)
    pub fn expected_count(&self, encoded: u64, k: usize, genome_length: u64) -> f64 {
        let n_positions = (genome_length as f64) - (k as f64) + 1.0;
        if n_positions <= 0.0 {
            return 0.0;
        }
        let mut prob = 1.0;
        let mut val = encoded;
        for _ in 0..k {
            let base = (val & 0b11) as usize;
            prob *= self.freqs[base];
            val >>= 2;
        }
        n_positions * prob
    }

    /// Significance threshold: λ + z × sqrt(max(λ, 1.0))
    pub fn significance_threshold(&self, encoded: u64, k: usize, genome_length: u64, z_score: f64) -> f64 {
        let lambda = self.expected_count(encoded, k, genome_length);
        lambda + z_score * lambda.max(1.0).sqrt()
    }

    /// Expected count for canonical k-mer: sum of forward + RC expected counts.
    /// For palindromes (forward == RC), returns just the forward expected count.
    pub fn expected_count_canonical(&self, encoded: u64, k: usize, genome_length: u64) -> f64 {
        let kmer = KmerU64 { encoded, len: k as u8 };
        let rc = kmer.reverse_complement();
        let forward_lambda = self.expected_count(encoded, k, genome_length);
        if rc.encoded == encoded {
            // Palindrome — don't double-count
            forward_lambda
        } else {
            forward_lambda + self.expected_count(rc.encoded, k, genome_length)
        }
    }

    /// Significance threshold for canonical k-mer.
    pub fn significance_threshold_canonical(&self, encoded: u64, k: usize, genome_length: u64, z_score: f64) -> f64 {
        let lambda = self.expected_count_canonical(encoded, k, genome_length);
        lambda + z_score * lambda.max(1.0).sqrt()
    }
}

impl Vocabulary {
    pub fn new() -> Self {
        Vocabulary {
            entries: FxHashMap::default(),
        }
    }

    pub fn insert(&mut self, kmer: KmerU64, count: u32) {
        self.entries
            .insert((kmer.encoded, kmer.len), VocabEntry { kmer, count });
    }

    #[inline]
    pub fn contains(&self, kmer: &KmerU64) -> bool {
        self.entries.contains_key(&(kmer.encoded, kmer.len))
    }

    pub fn get(&self, kmer: &KmerU64) -> Option<&VocabEntry> {
        self.entries.get(&(kmer.encoded, kmer.len))
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Get all k-mers of a specific length.
    pub fn kmers_of_length(&self, k: usize) -> Vec<&VocabEntry> {
        self.entries
            .values()
            .filter(|e| e.kmer.len as usize == k)
            .collect()
    }
}

/// Determine adaptive maximum k based on genome length.
pub fn find_max_k(genome_length: u64) -> usize {
    if genome_length == 0 {
        return 1;
    }
    let log4 = (genome_length as f64).log2() / 2.0;
    let max_k = (log4 + 2.0).floor() as usize;
    max_k.clamp(4, 32)
}

/// Count k-mers of a given length across all regions, in parallel.
/// Each region is counted independently on its own thread, then merged.
fn count_kmers_at_k(reader: &FastaReader, k: usize) -> FxHashMap<u64, u32> {
    if k > 32 {
        return FxHashMap::default();
    }

    let per_region: Vec<FxHashMap<u64, u32>> = reader
        .regions()
        .par_iter()
        .map(|region| count_kmers_raw(reader.region_bytes(region), k))
        .collect();

    let mut merged = FxHashMap::default();
    for counts in per_region {
        for (kmer, count) in counts {
            *merged.entry(kmer).or_insert(0) += count;
        }
    }

    merged
}

/// Count k-mers in raw FASTA bytes (with embedded newlines).
/// Newlines are skipped without resetting the sliding window.
/// Non-ACGT bases (N, etc.) reset the window.
fn count_kmers_raw(data: &[u8], k: usize) -> FxHashMap<u64, u32> {
    let mut counts = FxHashMap::default();
    if data.len() < k {
        return counts;
    }

    let mask = if k >= 32 {
        u64::MAX
    } else {
        (1u64 << (k * 2)) - 1
    };

    let mut encoded: u64 = 0;
    let mut valid_bases: usize = 0;

    for &byte in data {
        // Skip whitespace — don't reset the window
        if byte == b'\n' || byte == b'\r' {
            continue;
        }
        if let Some(bits) = encode::encode_base(byte) {
            encoded = ((encoded << 2) | bits as u64) & mask;
            valid_bases += 1;
        } else {
            // N or other invalid base — reset
            valid_bases = 0;
            encoded = 0;
            continue;
        }
        if valid_bases >= k {
            *counts.entry(encoded).or_insert(0) += 1;
        }
    }

    counts
}

/// Merge k-mer counts so that a k-mer and its reverse complement share one entry.
fn canonicalize_counts(counts: FxHashMap<u64, u32>, k: usize) -> FxHashMap<u64, u32> {
    let mut canonical: FxHashMap<u64, u32> = FxHashMap::default();
    for (&encoded, &count) in &counts {
        let kmer = KmerU64 { encoded, len: k as u8 };
        let canon = kmer.canonical();
        *canonical.entry(canon.encoded).or_insert(0) += count;
    }
    canonical
}

/// Build the full vocabulary across all k values from 1 to max_k.
/// Uses composition-aware Poisson significance testing with the given z_score.
pub fn build_vocabulary(reader: &FastaReader, max_k: usize, canonical: bool, z_score: f64) -> Vocabulary {
    let genome_length = reader.genome_length();
    let mut vocab = Vocabulary::new();

    // k=1: count once (always included, also used to compute base frequencies)
    let counts_1 = count_kmers_at_k(reader, 1);

    // Compute base frequencies from raw k=1 counts
    let base_freqs = BaseFreqs::from_counts(&counts_1);
    tracing::info!(
        "Base frequencies: A={:.3} C={:.3} G={:.3} T={:.3}",
        base_freqs.freqs[0],
        base_freqs.freqs[1],
        base_freqs.freqs[2],
        base_freqs.freqs[3],
    );

    if canonical {
        // Merge complement pairs: A+T -> A, C+G -> C
        let a = KmerU64::from_slice(b"A").unwrap();
        let t = KmerU64::from_slice(b"T").unwrap();
        let c = KmerU64::from_slice(b"C").unwrap();
        let g = KmerU64::from_slice(b"G").unwrap();
        let a_count = counts_1.get(&a.encoded).copied().unwrap_or(0)
            + counts_1.get(&t.encoded).copied().unwrap_or(0);
        let c_count = counts_1.get(&c.encoded).copied().unwrap_or(0)
            + counts_1.get(&g.encoded).copied().unwrap_or(0);
        vocab.insert(a, a_count);
        vocab.insert(c, c_count);
    } else {
        for &base in &[b'A', b'C', b'G', b'T'] {
            let kmer = KmerU64::from_slice(&[base]).unwrap();
            let count = counts_1.get(&kmer.encoded).copied().unwrap_or(0);
            vocab.insert(kmer, count);
        }
    }

    // k=2 to max_k: count and filter using per-k-mer significance threshold
    for k in 2..=max_k {
        let counts = count_kmers_at_k(reader, k);

        let filtered_counts = if canonical {
            canonicalize_counts(counts, k)
        } else {
            counts
        };

        let mut found = 0;
        let mut min_thresh = f64::MAX;
        let mut max_thresh = f64::MIN;

        for (&encoded, &count) in &filtered_counts {
            let threshold = if canonical {
                base_freqs.significance_threshold_canonical(encoded, k, genome_length, z_score)
            } else {
                base_freqs.significance_threshold(encoded, k, genome_length, z_score)
            };

            if threshold < min_thresh {
                min_thresh = threshold;
            }
            if threshold > max_thresh {
                max_thresh = threshold;
            }

            if (count as f64) > threshold {
                let kmer = KmerU64 {
                    encoded,
                    len: k as u8,
                };
                vocab.insert(kmer, count);
                found += 1;
            }
        }

        tracing::debug!(
            "k={}: {} passed filter (threshold range {:.1}..{:.1}, z={}, distinct={})",
            k,
            found,
            if min_thresh == f64::MAX { 0.0 } else { min_thresh },
            if max_thresh == f64::MIN { 0.0 } else { max_thresh },
            z_score,
            filtered_counts.len()
        );

        if found == 0 {
            tracing::info!("No k-mers passed filter at k={}, stopping", k);
            break;
        }
    }

    vocab
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_max_k() {
        assert!(find_max_k(1000) >= 4);
        assert!(find_max_k(3_000_000_000) <= 32);
        assert_eq!(find_max_k(0), 1);
    }

    #[test]
    fn test_count_kmers_raw() {
        // Simulate raw FASTA bytes (with newline mid-sequence)
        let data = b"ACGTAC\nGT";
        let counts = count_kmers_raw(data, 2);

        let ac = KmerU64::from_slice(b"AC").unwrap();
        let cg = KmerU64::from_slice(b"CG").unwrap();
        let gt = KmerU64::from_slice(b"GT").unwrap();
        let ta = KmerU64::from_slice(b"TA").unwrap();

        // ACGTACGT -> AC CG GT TA AC CG GT
        assert_eq!(counts[&ac.encoded], 2);
        assert_eq!(counts[&cg.encoded], 2);
        assert_eq!(counts[&gt.encoded], 2);
        assert_eq!(counts[&ta.encoded], 1);
    }

    #[test]
    fn test_count_kmers_raw_with_n() {
        // N should reset the window
        let data = b"ACNGT";
        let counts = count_kmers_raw(data, 2);

        let ac = KmerU64::from_slice(b"AC").unwrap();
        let gt = KmerU64::from_slice(b"GT").unwrap();

        assert_eq!(counts.get(&ac.encoded).copied().unwrap_or(0), 1);
        assert_eq!(counts.get(&gt.encoded).copied().unwrap_or(0), 1);
        // "CN" and "NG" should NOT exist
        assert_eq!(counts.len(), 2);
    }

    #[test]
    fn test_vocabulary_lookup() {
        let mut vocab = Vocabulary::new();
        let kmer = KmerU64::from_slice(b"ACG").unwrap();
        vocab.insert(kmer, 10);
        assert!(vocab.contains(&kmer));
        assert_eq!(vocab.get(&kmer).unwrap().count, 10);

        let other = KmerU64::from_slice(b"TTT").unwrap();
        assert!(!vocab.contains(&other));
    }

    #[test]
    fn test_canonical_vocabulary() {
        // Test canonicalize_counts directly
        let mut counts = FxHashMap::default();
        let ac = KmerU64::from_slice(b"AC").unwrap();
        let gt = KmerU64::from_slice(b"GT").unwrap();
        counts.insert(ac.encoded, 5);
        counts.insert(gt.encoded, 3);

        let canon = canonicalize_counts(counts, 2);
        // AC and GT are RC pairs, should merge
        let ac_canon = ac.canonical();
        assert_eq!(canon.len(), 1); // merged into one
        assert_eq!(canon[&ac_canon.encoded], 8); // 5 + 3
    }

    #[test]
    fn test_base_freqs_uniform() {
        let mut counts = FxHashMap::default();
        counts.insert(0, 100); // A
        counts.insert(1, 100); // C
        counts.insert(2, 100); // G
        counts.insert(3, 100); // T
        let bf = BaseFreqs::from_counts(&counts);
        for &f in &bf.freqs {
            assert!((f - 0.25).abs() < 1e-10);
        }
    }

    #[test]
    fn test_base_freqs_biased() {
        // Simulate AT-rich genome: A=42%, T=41%, C=8.5%, G=8.5%
        let mut counts = FxHashMap::default();
        counts.insert(0, 4200); // A
        counts.insert(1, 850);  // C
        counts.insert(2, 850);  // G
        counts.insert(3, 4100); // T
        let bf = BaseFreqs::from_counts(&counts);
        assert!((bf.freqs[0] - 0.42).abs() < 0.001);
        assert!((bf.freqs[1] - 0.085).abs() < 0.001);
        assert!((bf.freqs[2] - 0.085).abs() < 0.001);
        assert!((bf.freqs[3] - 0.41).abs() < 0.001);
    }

    #[test]
    fn test_expected_count_biased() {
        let mut counts = FxHashMap::default();
        counts.insert(0, 4200); // A
        counts.insert(1, 850);  // C
        counts.insert(2, 850);  // G
        counts.insert(3, 4100); // T
        let bf = BaseFreqs::from_counts(&counts);

        let genome_length = 100_000u64;

        // AAAA: expected = (100000-3) * 0.42^4
        let aaaa = KmerU64::from_slice(b"AAAA").unwrap();
        let expected = bf.expected_count(aaaa.encoded, 4, genome_length);
        let manual = (genome_length as f64 - 3.0) * 0.42_f64.powi(4);
        assert!((expected - manual).abs() / manual < 0.01);

        // CCCC: expected = (100000-3) * 0.085^4
        let cccc = KmerU64::from_slice(b"CCCC").unwrap();
        let expected_cccc = bf.expected_count(cccc.encoded, 4, genome_length);
        let manual_cccc = (genome_length as f64 - 3.0) * 0.085_f64.powi(4);
        assert!((expected_cccc - manual_cccc).abs() / manual_cccc < 0.01);

        // AAAA should have much higher expected count than CCCC
        assert!(expected > expected_cccc * 10.0);
    }

    #[test]
    fn test_significance_threshold_floor() {
        // When lambda is very small, max(lambda, 1.0) floor ensures minimum threshold
        let mut counts = FxHashMap::default();
        counts.insert(0, 4200);
        counts.insert(1, 850);
        counts.insert(2, 850);
        counts.insert(3, 4100);
        let bf = BaseFreqs::from_counts(&counts);

        // A rare k-mer at high k should get a floor from max(lambda, 1.0)
        let ccccggggg = KmerU64::from_slice(b"CCCCGGGGG").unwrap();
        let threshold = bf.significance_threshold(ccccggggg.encoded, 9, 100_000, 3.0);
        // With z=3.0 and lambda near 0, threshold should be at least 3.0 (z * sqrt(1.0))
        assert!(threshold >= 3.0);
    }

    #[test]
    fn test_significance_threshold_high_lambda() {
        // When lambda is large, threshold scales proportionally
        let mut counts = FxHashMap::default();
        counts.insert(0, 2500);
        counts.insert(1, 2500);
        counts.insert(2, 2500);
        counts.insert(3, 2500);
        let bf = BaseFreqs::from_counts(&counts);

        // AAAA in a uniform 1M genome: lambda = (1000000-3) * 0.25^4 ≈ 3906
        let aaaa = KmerU64::from_slice(b"AAAA").unwrap();
        let lambda = bf.expected_count(aaaa.encoded, 4, 1_000_000);
        let threshold = bf.significance_threshold(aaaa.encoded, 4, 1_000_000, 3.0);

        // threshold = lambda + 3 * sqrt(lambda) ≈ 3906 + 3*62.5 ≈ 4094
        assert!(threshold > lambda);
        assert!((threshold - lambda - 3.0 * lambda.sqrt()).abs() < 1.0);
    }
}
