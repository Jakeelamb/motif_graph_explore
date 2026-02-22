/// Multi-genome comparative analysis.
///
/// Computes per-genome statistics from assembly index results and produces
/// side-by-side comparison output in CSV and JSON formats.

use crate::assembly::AssemblyResult;
use crate::graph::DagStats;
use crate::output::Metadata;
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::io::Write;
use std::path::Path;

/// Statistics for a single genome.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GenomeStats {
    pub label: String,
    pub genome_length: u64,
    pub vocab_size: usize,
    pub max_assembly_index: u32,
    pub mean_assembly_index: f64,
    pub median_assembly_index: u32,
    pub node_count: usize,
    pub edge_count: usize,
    pub max_depth: usize,
    pub ai_distribution: Vec<(u32, usize)>,
    pub level_distribution: Vec<(usize, usize)>,
}

/// Result of comparing multiple genomes.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ComparisonResult {
    pub genomes: Vec<GenomeStats>,
    pub genome_count: usize,
}

/// Compute statistics for a single genome from its assembly results.
pub fn compute_genome_stats(
    label: &str,
    result: &AssemblyResult,
    dag_stats: &DagStats,
    metadata: &Metadata,
) -> GenomeStats {
    let mut ai_values: Vec<u32> = result
        .nodes
        .iter()
        .map(|n| n.assembly_index)
        .collect();
    ai_values.sort_unstable();

    let mean_ai = if ai_values.is_empty() {
        0.0
    } else {
        ai_values.iter().map(|&v| v as f64).sum::<f64>() / ai_values.len() as f64
    };

    let median_ai = if ai_values.is_empty() {
        0
    } else {
        ai_values[ai_values.len() / 2]
    };

    // AI distribution: count of nodes at each AI value
    let mut ai_dist: Vec<(u32, usize)> = Vec::new();
    if !ai_values.is_empty() {
        let mut current_ai = ai_values[0];
        let mut count = 0usize;
        for &ai in &ai_values {
            if ai == current_ai {
                count += 1;
            } else {
                ai_dist.push((current_ai, count));
                current_ai = ai;
                count = 1;
            }
        }
        ai_dist.push((current_ai, count));
    }

    GenomeStats {
        label: label.to_string(),
        genome_length: metadata.genome_length,
        vocab_size: metadata.vocab_size,
        max_assembly_index: metadata.max_assembly_index,
        mean_assembly_index: mean_ai,
        median_assembly_index: median_ai,
        node_count: dag_stats.node_count,
        edge_count: dag_stats.edge_count,
        max_depth: dag_stats.max_depth,
        ai_distribution: ai_dist,
        level_distribution: dag_stats.nodes_per_level.clone(),
    }
}

/// Compare multiple genomes by collecting their stats.
pub fn compare_genomes(entries: Vec<GenomeStats>) -> ComparisonResult {
    let genome_count = entries.len();
    ComparisonResult {
        genomes: entries,
        genome_count,
    }
}

/// Write comparison results as CSV.
pub fn write_comparison_csv(comparison: &ComparisonResult, path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(path)
        .with_context(|| format!("Failed to create comparison CSV: {}", path.display()))?;

    writeln!(
        f,
        "label,genome_length,vocab_size,max_ai,mean_ai,median_ai,node_count,edge_count,max_depth"
    )?;

    for g in &comparison.genomes {
        writeln!(
            f,
            "{},{},{},{},{:.4},{},{},{},{}",
            g.label,
            g.genome_length,
            g.vocab_size,
            g.max_assembly_index,
            g.mean_assembly_index,
            g.median_assembly_index,
            g.node_count,
            g.edge_count,
            g.max_depth,
        )?;
    }

    Ok(())
}

/// Write comparison results as JSON.
pub fn write_comparison_json(comparison: &ComparisonResult, path: &Path) -> Result<()> {
    let json = serde_json::to_string_pretty(comparison)?;
    std::fs::write(path, json)
        .with_context(|| format!("Failed to write comparison JSON: {}", path.display()))?;
    Ok(())
}

/// Print a summary comparison table to stdout.
pub fn print_comparison_table(comparison: &ComparisonResult) {
    println!(
        "{:<30} {:>12} {:>8} {:>6} {:>8} {:>6} {:>8} {:>8} {:>5}",
        "Genome", "Length", "Vocab", "MaxAI", "MeanAI", "MedAI", "Nodes", "Edges", "Depth"
    );
    println!("{}", "-".repeat(105));

    for g in &comparison.genomes {
        println!(
            "{:<30} {:>12} {:>8} {:>6} {:>8.2} {:>6} {:>8} {:>8} {:>5}",
            g.label,
            g.genome_length,
            g.vocab_size,
            g.max_assembly_index,
            g.mean_assembly_index,
            g.median_assembly_index,
            g.node_count,
            g.edge_count,
            g.max_depth,
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::{AssemblyNode, AssemblyResult};
    use crate::encode::KmerU64;
    use crate::graph::DagStats;
    use crate::output::Metadata;

    fn make_test_result() -> (AssemblyResult, DagStats, Metadata) {
        let nodes = vec![
            AssemblyNode {
                kmer: KmerU64::from_slice(b"A").unwrap(),
                assembly_index: 0,
                count: 100,
                split_pos: None,
            },
            AssemblyNode {
                kmer: KmerU64::from_slice(b"C").unwrap(),
                assembly_index: 0,
                count: 100,
                split_pos: None,
            },
            AssemblyNode {
                kmer: KmerU64::from_slice(b"AC").unwrap(),
                assembly_index: 3,
                count: 50,
                split_pos: Some(1),
            },
        ];
        let result = AssemblyResult { nodes };

        let dag_stats = DagStats {
            node_count: 3,
            edge_count: 2,
            max_depth: 2,
            nodes_per_level: vec![(1, 2), (2, 1)],
        };

        let metadata = Metadata {
            genome_file: "test.fasta".to_string(),
            genome_length: 1000,
            max_k: 4,
            vocab_size: 3,
            max_assembly_index: 3,
            node_count: 3,
            edge_count: 2,
            elapsed_secs: 0.1,
            canonical: false,
            z_score: 3.0,
        };

        (result, dag_stats, metadata)
    }

    #[test]
    fn test_compute_genome_stats() {
        let (result, dag_stats, metadata) = make_test_result();
        let stats = compute_genome_stats("test_genome", &result, &dag_stats, &metadata);

        assert_eq!(stats.label, "test_genome");
        assert_eq!(stats.genome_length, 1000);
        assert_eq!(stats.max_assembly_index, 3);
        assert_eq!(stats.median_assembly_index, 0);
        assert!((stats.mean_assembly_index - 1.0).abs() < 0.01); // (0+0+3)/3 = 1.0
        assert_eq!(stats.node_count, 3);
        assert_eq!(stats.edge_count, 2);

        // AI distribution: two 0s, one 3
        assert_eq!(stats.ai_distribution, vec![(0, 2), (3, 1)]);
    }

    #[test]
    fn test_write_comparison_csv() {
        let (result, dag_stats, metadata) = make_test_result();
        let stats = compute_genome_stats("genome1", &result, &dag_stats, &metadata);
        let comparison = compare_genomes(vec![stats]);

        let tmp = std::env::temp_dir().join("test_comparison.csv");
        write_comparison_csv(&comparison, &tmp).unwrap();

        let content = std::fs::read_to_string(&tmp).unwrap();
        assert!(content.contains("label,genome_length,vocab_size"));
        assert!(content.contains("genome1,1000,3,3"));
        std::fs::remove_file(tmp).ok();
    }

    #[test]
    fn test_compare_genomes() {
        let (result, dag_stats, metadata) = make_test_result();
        let stats1 = compute_genome_stats("genome1", &result, &dag_stats, &metadata);
        let stats2 = compute_genome_stats("genome2", &result, &dag_stats, &metadata);
        let comparison = compare_genomes(vec![stats1, stats2]);

        assert_eq!(comparison.genome_count, 2);
        assert_eq!(comparison.genomes.len(), 2);
        assert_eq!(comparison.genomes[0].label, "genome1");
        assert_eq!(comparison.genomes[1].label, "genome2");
    }
}
