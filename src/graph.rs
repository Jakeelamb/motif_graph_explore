/// Assembly DAG construction from assembly index results.
///
/// Nodes are vocabulary k-mers, edges represent optimal decompositions
/// (parent -> left child, parent -> right child).

use crate::assembly::AssemblyResult;
use petgraph::graph::{DiGraph, NodeIndex};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::fmt;

/// Which half of a split an edge represents.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Side {
    Left,
    Right,
}

impl fmt::Display for Side {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Side::Left => write!(f, "left"),
            Side::Right => write!(f, "right"),
        }
    }
}

/// Attributes stored on each graph node.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NodeAttr {
    pub sequence: String,
    pub assembly_index: u32,
    pub frequency: u32,
    pub level: u8,
}

/// Attributes stored on each graph edge.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct EdgeAttr {
    pub side: Side,
    pub split_pos: usize,
}

/// The assembly DAG.
pub type AssemblyDag = DiGraph<NodeAttr, EdgeAttr>;

/// Statistics about the assembly DAG.
#[derive(Debug, Serialize, Deserialize)]
pub struct DagStats {
    pub node_count: usize,
    pub edge_count: usize,
    pub max_depth: usize,
    pub nodes_per_level: Vec<(usize, usize)>,
}

/// Build the assembly DAG from computation results.
pub fn build_dag(result: &AssemblyResult) -> AssemblyDag {
    let mut graph = DiGraph::new();
    let mut node_map: FxHashMap<(u64, u8), NodeIndex> = FxHashMap::default();

    // Create all nodes
    for node in &result.nodes {
        let idx = graph.add_node(NodeAttr {
            sequence: node.kmer.to_string(),
            assembly_index: node.assembly_index,
            frequency: node.count,
            level: node.kmer.len,
        });
        node_map.insert((node.kmer.encoded, node.kmer.len), idx);
    }

    // Add edges for optimal splits
    for node in &result.nodes {
        if let Some(split_pos) = node.split_pos {
            let parent_idx = node_map[&(node.kmer.encoded, node.kmer.len)];
            let (left, right) = node.kmer.split_at(split_pos);

            if let Some(&left_idx) = node_map.get(&(left.encoded, left.len)) {
                graph.add_edge(parent_idx, left_idx, EdgeAttr { side: Side::Left, split_pos });
            }

            if let Some(&right_idx) = node_map.get(&(right.encoded, right.len)) {
                graph.add_edge(parent_idx, right_idx, EdgeAttr { side: Side::Right, split_pos });
            }
        }
    }

    graph
}

/// Compute statistics about the assembly DAG.
pub fn dag_stats(dag: &AssemblyDag) -> DagStats {
    let mut level_counts: FxHashMap<usize, usize> = FxHashMap::default();
    let mut max_level = 0usize;

    for attr in dag.node_weights() {
        let level = attr.level as usize;
        *level_counts.entry(level).or_insert(0) += 1;
        max_level = max_level.max(level);
    }

    let mut nodes_per_level: Vec<(usize, usize)> = level_counts.into_iter().collect();
    nodes_per_level.sort_by_key(|&(level, _)| level);

    DagStats {
        node_count: dag.node_count(),
        edge_count: dag.edge_count(),
        max_depth: max_level,
        nodes_per_level,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::compute_assembly_index;
    use crate::encode::KmerU64;
    use crate::kmer::Vocabulary;

    #[test]
    fn test_build_dag() {
        let mut vocab = Vocabulary::new();
        for &base in &[b'A', b'C', b'G', b'T'] {
            vocab.insert(KmerU64::from_slice(&[base]).unwrap(), 100);
        }
        vocab.insert(KmerU64::from_slice(b"AC").unwrap(), 50);

        let result = compute_assembly_index(&vocab, 2, false);
        let dag = build_dag(&result);

        assert_eq!(dag.node_count(), 5);
        assert_eq!(dag.edge_count(), 2);
    }

    #[test]
    fn test_dag_is_acyclic() {
        let mut vocab = Vocabulary::new();
        for &base in &[b'A', b'C', b'G', b'T'] {
            vocab.insert(KmerU64::from_slice(&[base]).unwrap(), 100);
        }
        vocab.insert(KmerU64::from_slice(b"AC").unwrap(), 50);
        vocab.insert(KmerU64::from_slice(b"CG").unwrap(), 50);
        vocab.insert(KmerU64::from_slice(b"ACG").unwrap(), 25);

        let result = compute_assembly_index(&vocab, 3, false);
        let dag = build_dag(&result);

        assert!(petgraph::algo::toposort(&dag, None).is_ok());
    }

    #[test]
    fn test_dag_stats() {
        let mut vocab = Vocabulary::new();
        for &base in &[b'A', b'C', b'G', b'T'] {
            vocab.insert(KmerU64::from_slice(&[base]).unwrap(), 100);
        }
        vocab.insert(KmerU64::from_slice(b"AC").unwrap(), 50);

        let result = compute_assembly_index(&vocab, 2, false);
        let dag = build_dag(&result);
        let stats = dag_stats(&dag);

        assert_eq!(stats.node_count, 5);
        assert_eq!(stats.max_depth, 2);
    }
}
