# Assembly Theory for Genomic Sequences

## What This Computes

Assembly Index (AI) measures the minimum number of joining operations needed to construct
a string from its basic alphabet, reusing previously constructed substrings. For DNA:

- Base alphabet: {A, T, C, G} — AI = 0
- Dinucleotides (AT, CG, ...): one join → AI = 1
- Longer k-mers: recursively built from shorter pieces already in the vocabulary

The key insight: substrings that appear more often than random expectation form a
"vocabulary" that the genome reuses. The assembly DAG captures how complex motifs
are built from simpler ones through concatenation.

## Why It Matters

Assembly Theory (Cronin et al.) provides a complexity measure grounded in constructive
steps rather than information entropy. Applied to genomes, it reveals:

- Which sequence motifs are compositionally complex
- How the genome's repetitive vocabulary is hierarchically structured
- Shared construction pathways between different motifs

## Algorithm

1. Extract k-mer vocabulary: sequences occurring above random expectation
2. Bottom-up dynamic programming: for each k-mer, find the split into two
   vocabulary members that minimizes the total pathway (set of unique
   intermediates needed)
3. Build the assembly DAG: nodes are vocabulary k-mers, edges are optimal splits
4. Export for exploration

## Design Choices

- **Rust** for performance on genome-scale data
- **2-bit encoding** packs k-mers (k≤32) into u64 for O(1) hashing
- **Level-parallel** DP: all k-mers of the same length are independent
- **Pathway tracking** via sets, not just counts, to handle shared intermediates
