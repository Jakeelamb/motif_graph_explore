/// Assembly Index computation via bottom-up dynamic programming.
///
/// For each k-mer in the vocabulary, finds the minimum-cost way to construct it
/// by concatenating two shorter vocabulary members. The "cost" is the size of the
/// pathway — the set of all unique intermediate k-mers needed in the construction.
///
/// Base alphabet {A, C, G, T} has AI = 0. Each level builds on the previous.
///
/// Pathway sets are tracked during DP computation (needed for correct union
/// operations) but not stored in the final output — only AI and split position
/// are kept per node.

use crate::encode::KmerU64;
use crate::kmer::Vocabulary;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};

/// Result of assembly index computation for a single k-mer.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AssemblyNode {
    pub kmer: KmerU64,
    pub assembly_index: u32,
    pub count: u32,
    /// Optimal split position (None for base alphabet).
    pub split_pos: Option<usize>,
}

/// Full result of the assembly index computation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AssemblyResult {
    pub nodes: Vec<AssemblyNode>,
}

impl AssemblyResult {
    pub fn max_assembly_index(&self) -> u32 {
        self.nodes
            .iter()
            .map(|n| n.assembly_index)
            .max()
            .unwrap_or(0)
    }

    pub fn get(&self, kmer: &KmerU64) -> Option<&AssemblyNode> {
        self.nodes.iter().find(|n| n.kmer == *kmer)
    }
}

/// Compute assembly index for all k-mers in the vocabulary.
pub fn compute_assembly_index(vocabulary: &Vocabulary, max_k: usize, canonical: bool) -> AssemblyResult {
    // Pathway storage for the DP. Only lives during computation.
    // Key: (encoded, len), Value: (ai, split_pos, pathway_kmers)
    let mut computed: FxHashMap<(u64, u8), (u32, Option<usize>, Vec<KmerU64>)> =
        FxHashMap::default();

    let mut result_nodes = Vec::new();

    // Level 0: base alphabet — AI = 0
    let base_alphabet: &[u8] = if canonical {
        &[b'A', b'C']
    } else {
        &[b'A', b'C', b'G', b'T']
    };
    for &base in base_alphabet {
        let kmer = KmerU64::from_slice(&[base]).unwrap();
        if vocabulary.contains(&kmer) {
            let count = vocabulary.get(&kmer).map(|e| e.count).unwrap_or(0);
            result_nodes.push(AssemblyNode {
                kmer,
                assembly_index: 0,
                count,
                split_pos: None,
            });
            computed.insert((kmer.encoded, kmer.len), (0, None, vec![kmer]));
        }
    }

    // Process levels bottom-up: k=2, k=3, ..., k=max_k
    for k in 2..=max_k {
        let level_kmers = vocabulary.kmers_of_length(k);
        if level_kmers.is_empty() {
            continue;
        }

        // All k-mers at the same level are independent — parallel.
        // Only reads entries with len < k, which are all finalized.
        let computed_ref = &computed;

        let results: Vec<(AssemblyNode, Vec<KmerU64>)> = level_kmers
            .par_iter()
            .filter_map(|entry| {
                let kmer = entry.kmer;
                let mut best_ai = u32::MAX;
                let mut best_split: Option<usize> = None;

                // Reuse one HashSet across all split positions for this k-mer
                let mut pathway_set: FxHashSet<(u64, u8)> = FxHashSet::default();

                for p in 1..k {
                    let (left, right) = kmer.split_at(p);

                    // When canonical, look up canonical forms of split halves
                    let (left_lookup, right_lookup) = if canonical {
                        (left.canonical(), right.canonical())
                    } else {
                        (left, right)
                    };

                    if !vocabulary.contains(&left_lookup) || !vocabulary.contains(&right_lookup) {
                        continue;
                    }

                    let left_data = match computed_ref.get(&(left_lookup.encoded, left_lookup.len)) {
                        Some(d) => d,
                        None => continue,
                    };
                    let right_data = match computed_ref.get(&(right_lookup.encoded, right_lookup.len)) {
                        Some(d) => d,
                        None => continue,
                    };

                    pathway_set.clear();
                    for km in &left_data.2 {
                        pathway_set.insert((km.encoded, km.len));
                    }
                    for km in &right_data.2 {
                        pathway_set.insert((km.encoded, km.len));
                    }
                    pathway_set.insert((kmer.encoded, kmer.len));

                    let ai = pathway_set.len() as u32;
                    if ai < best_ai {
                        best_ai = ai;
                        best_split = Some(p);
                    }
                }

                // Rebuild the best pathway once (for storage in the DP table)
                let split = best_split?;
                let (left, right) = kmer.split_at(split);
                let (left_lookup, right_lookup) = if canonical {
                    (left.canonical(), right.canonical())
                } else {
                    (left, right)
                };
                let left_data = computed_ref.get(&(left_lookup.encoded, left_lookup.len)).unwrap();
                let right_data = computed_ref.get(&(right_lookup.encoded, right_lookup.len)).unwrap();

                pathway_set.clear();
                for km in &left_data.2 {
                    pathway_set.insert((km.encoded, km.len));
                }
                for km in &right_data.2 {
                    pathway_set.insert((km.encoded, km.len));
                }
                pathway_set.insert((kmer.encoded, kmer.len));

                let pathway: Vec<KmerU64> = pathway_set
                    .iter()
                    .map(|&(enc, len)| KmerU64 {
                        encoded: enc,
                        len,
                    })
                    .collect();

                Some((
                    AssemblyNode {
                        kmer,
                        assembly_index: best_ai,
                        count: entry.count,
                        split_pos: Some(split),
                    },
                    pathway,
                ))
            })
            .collect();

        // Merge: move pathways into computed (no clone), push nodes into results
        for (node, pathway) in results {
            computed.insert(
                (node.kmer.encoded, node.kmer.len),
                (node.assembly_index, node.split_pos, pathway),
            );
            result_nodes.push(node);
        }

        tracing::debug!(
            "Level k={}: {} assembly indices computed",
            k,
            level_kmers.len()
        );
    }

    AssemblyResult {
        nodes: result_nodes,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::Vocabulary;

    fn make_test_vocabulary() -> Vocabulary {
        let mut vocab = Vocabulary::new();
        for &base in &[b'A', b'C', b'G', b'T'] {
            vocab.insert(KmerU64::from_slice(&[base]).unwrap(), 100);
        }
        vocab.insert(KmerU64::from_slice(b"AC").unwrap(), 50);
        vocab.insert(KmerU64::from_slice(b"GT").unwrap(), 50);
        vocab.insert(KmerU64::from_slice(b"CG").unwrap(), 50);
        vocab.insert(KmerU64::from_slice(b"ACG").unwrap(), 25);
        vocab
    }

    #[test]
    fn test_base_alphabet_ai_zero() {
        let vocab = make_test_vocabulary();
        let result = compute_assembly_index(&vocab, 3, false);

        let a = result.get(&KmerU64::from_slice(b"A").unwrap()).unwrap();
        assert_eq!(a.assembly_index, 0);
        assert!(a.split_pos.is_none());
    }

    #[test]
    fn test_dinucleotide_ai() {
        let vocab = make_test_vocabulary();
        let result = compute_assembly_index(&vocab, 3, false);

        let ac = result.get(&KmerU64::from_slice(b"AC").unwrap()).unwrap();
        // pathway = {A} ∪ {C} ∪ {AC} = {A, C, AC} -> AI = 3
        assert_eq!(ac.assembly_index, 3);
        assert!(ac.split_pos.is_some());
    }

    #[test]
    fn test_trinucleotide_ai() {
        let vocab = make_test_vocabulary();
        let result = compute_assembly_index(&vocab, 3, false);

        let acg = result.get(&KmerU64::from_slice(b"ACG").unwrap()).unwrap();
        // A+CG: {A} ∪ {C,G,CG} ∪ {ACG} = {A,C,G,CG,ACG} -> AI=5
        // AC+G: {A,C,AC} ∪ {G} ∪ {ACG} = {A,C,AC,G,ACG} -> AI=5
        assert_eq!(acg.assembly_index, 5);
    }

    #[test]
    fn test_no_valid_split() {
        let mut vocab = Vocabulary::new();
        for &base in &[b'A', b'C', b'G', b'T'] {
            vocab.insert(KmerU64::from_slice(&[base]).unwrap(), 100);
        }
        // ACG in vocab but no dinucleotides — can't split
        vocab.insert(KmerU64::from_slice(b"ACG").unwrap(), 25);

        let result = compute_assembly_index(&vocab, 3, false);
        assert!(result.get(&KmerU64::from_slice(b"ACG").unwrap()).is_none());
    }

    #[test]
    fn test_canonical_assembly_index() {
        let mut vocab = Vocabulary::new();
        // Canonical alphabet: only A and C
        vocab.insert(KmerU64::from_slice(b"A").unwrap(), 200);
        vocab.insert(KmerU64::from_slice(b"C").unwrap(), 200);
        // AC is canonical (AC < GT in encoding)
        vocab.insert(KmerU64::from_slice(b"AC").unwrap(), 100);

        let result = compute_assembly_index(&vocab, 2, true);

        let a = result.get(&KmerU64::from_slice(b"A").unwrap()).unwrap();
        assert_eq!(a.assembly_index, 0);

        let ac = result.get(&KmerU64::from_slice(b"AC").unwrap()).unwrap();
        // pathway = {A} ∪ {C} ∪ {AC} = 3
        assert_eq!(ac.assembly_index, 3);
    }
}
