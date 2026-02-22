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

/// Build the full vocabulary across all k values from 1 to max_k.
pub fn build_vocabulary(reader: &FastaReader, max_k: usize) -> Vocabulary {
    let genome_length = reader.genome_length();
    let mut vocab = Vocabulary::new();

    // k=1: count once, always include all 4 bases
    let counts_1 = count_kmers_at_k(reader, 1);
    for &base in &[b'A', b'C', b'G', b'T'] {
        let kmer = KmerU64::from_slice(&[base]).unwrap();
        let count = counts_1.get(&kmer.encoded).copied().unwrap_or(0);
        vocab.insert(kmer, count);
    }

    // k=2 to max_k: count and filter
    for k in 2..=max_k {
        let counts = count_kmers_at_k(reader, k);
        let threshold = random_expectation_count(genome_length, k);

        let mut found = 0;
        for (&encoded, &count) in &counts {
            if count as u64 > threshold {
                let kmer = KmerU64 {
                    encoded,
                    len: k as u8,
                };
                vocab.insert(kmer, count);
                found += 1;
            }
        }

        tracing::debug!(
            "k={}: {} passed filter (threshold={}, distinct={})",
            k,
            found,
            threshold,
            counts.len()
        );

        if found == 0 {
            tracing::info!("No k-mers passed filter at k={}, stopping", k);
            break;
        }
    }

    vocab
}

/// Expected count = genome_length / 4^k
fn random_expectation_count(genome_length: u64, k: usize) -> u64 {
    let four_to_k = 4u64.checked_pow(k as u32).unwrap_or(u64::MAX);
    genome_length / four_to_k
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
    fn test_random_expectation() {
        assert_eq!(random_expectation_count(1000, 2), 62);
        assert_eq!(random_expectation_count(1000, 5), 0);
    }
}
