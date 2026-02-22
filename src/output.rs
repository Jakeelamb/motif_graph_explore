/// Output serialization: GraphML, DOT, CSV, JSON, and binary formats.

use crate::assembly::AssemblyResult;
use crate::graph::{AssemblyDag, EdgeAttr, NodeAttr};
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;

/// Export format for the CLI.
#[derive(Clone, Debug, clap::ValueEnum)]
pub enum ExportFormat {
    Graphml,
    Dot,
    Csv,
    Json,
}

/// Run metadata.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Metadata {
    pub genome_file: String,
    pub genome_length: u64,
    pub max_k: usize,
    pub vocab_size: usize,
    pub max_assembly_index: u32,
    pub node_count: usize,
    pub edge_count: usize,
    pub elapsed_secs: f64,
    #[serde(default)]
    pub canonical: bool,
    #[serde(default = "default_z_score")]
    pub z_score: f64,
}

fn default_z_score() -> f64 {
    3.0
}

// --- DOT export (for graphviz: `dot -Tsvg file.dot -o file.svg`) ---

/// Write the assembly DAG as DOT (Graphviz).
pub fn write_dot(dag: &AssemblyDag, path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(path)
        .with_context(|| format!("Failed to create DOT file: {}", path.display()))?;

    writeln!(f, "digraph assembly_dag {{")?;
    writeln!(f, "    rankdir=BT;")?;
    writeln!(
        f,
        "    node [shape=box style=filled fontname=\"Courier\"];"
    )?;

    // Nodes
    for idx in dag.node_indices() {
        let attr = &dag[idx];
        writeln!(
            f,
            "    n{} [label=\"{}\" fillcolor=\"{}\" tooltip=\"AI={} freq={} k={}\"];",
            idx.index(),
            attr.sequence,
            ai_color(attr.assembly_index),
            attr.assembly_index,
            attr.frequency,
            attr.level,
        )?;
    }

    // Edges
    for eidx in dag.edge_indices() {
        let (s, t) = dag.edge_endpoints(eidx).unwrap();
        writeln!(f, "    n{} -> n{};", s.index(), t.index())?;
    }

    // Rank constraints: group nodes by level
    let mut by_level: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
    for idx in dag.node_indices() {
        by_level
            .entry(dag[idx].level)
            .or_default()
            .push(idx.index());
    }
    for (_, indices) in &by_level {
        write!(f, "    {{ rank=same;")?;
        for &idx in indices {
            write!(f, " n{};", idx)?;
        }
        writeln!(f, " }}")?;
    }

    writeln!(f, "}}")?;
    Ok(())
}

fn ai_color(ai: u32) -> &'static str {
    match ai {
        0 => "#66bb6a",
        1..=3 => "#42a5f5",
        4..=6 => "#ffa726",
        7..=9 => "#ef5350",
        _ => "#ab47bc",
    }
}

// --- GraphML export (for Gephi/Cytoscape) ---

/// Write the assembly DAG as GraphML.
pub fn write_graphml(dag: &AssemblyDag, path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(path)
        .with_context(|| format!("Failed to create GraphML file: {}", path.display()))?;

    writeln!(f, r#"<?xml version="1.0" encoding="UTF-8"?>"#)?;
    writeln!(
        f,
        r#"<graphml xmlns="http://graphml.graphstruct.org/graphml">"#
    )?;

    writeln!(
        f,
        r#"  <key id="sequence" for="node" attr.name="sequence" attr.type="string"/>"#
    )?;
    writeln!(
        f,
        r#"  <key id="assembly_index" for="node" attr.name="assembly_index" attr.type="int"/>"#
    )?;
    writeln!(
        f,
        r#"  <key id="frequency" for="node" attr.name="frequency" attr.type="int"/>"#
    )?;
    writeln!(
        f,
        r#"  <key id="level" for="node" attr.name="level" attr.type="int"/>"#
    )?;
    writeln!(
        f,
        r#"  <key id="side" for="edge" attr.name="side" attr.type="string"/>"#
    )?;
    writeln!(
        f,
        r#"  <key id="split_pos" for="edge" attr.name="split_pos" attr.type="int"/>"#
    )?;

    writeln!(f, r#"  <graph id="assembly_dag" edgedefault="directed">"#)?;

    for idx in dag.node_indices() {
        let attr = &dag[idx];
        writeln!(f, r#"    <node id="n{}">"#, idx.index())?;
        writeln!(f, r#"      <data key="sequence">{}</data>"#, attr.sequence)?;
        writeln!(
            f,
            r#"      <data key="assembly_index">{}</data>"#,
            attr.assembly_index
        )?;
        writeln!(f, r#"      <data key="frequency">{}</data>"#, attr.frequency)?;
        writeln!(f, r#"      <data key="level">{}</data>"#, attr.level)?;
        writeln!(f, r#"    </node>"#)?;
    }

    for eidx in dag.edge_indices() {
        let (source, target) = dag.edge_endpoints(eidx).unwrap();
        let attr = &dag[eidx];
        writeln!(
            f,
            r#"    <edge source="n{}" target="n{}">"#,
            source.index(),
            target.index()
        )?;
        writeln!(f, r#"      <data key="side">{}</data>"#, attr.side)?;
        writeln!(f, r#"      <data key="split_pos">{}</data>"#, attr.split_pos)?;
        writeln!(f, r#"    </edge>"#)?;
    }

    writeln!(f, r#"  </graph>"#)?;
    writeln!(f, r#"</graphml>"#)?;
    Ok(())
}

// --- CSV export ---

pub fn write_csv(result: &AssemblyResult, path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(path)
        .with_context(|| format!("Failed to create CSV file: {}", path.display()))?;

    writeln!(f, "kmer,length,assembly_index,frequency,split_pos")?;
    for node in &result.nodes {
        writeln!(
            f,
            "{},{},{},{},{}",
            node.kmer.to_string(),
            node.kmer.len,
            node.assembly_index,
            node.count,
            node.split_pos.map(|p| p.to_string()).unwrap_or_default()
        )?;
    }

    Ok(())
}

// --- JSON metadata ---

pub fn write_json(metadata: &Metadata, path: &Path) -> Result<()> {
    let json = serde_json::to_string_pretty(metadata)?;
    std::fs::write(path, json)
        .with_context(|| format!("Failed to write JSON: {}", path.display()))?;
    Ok(())
}

// --- Binary serialization (for export/info subcommands) ---

#[derive(Serialize, Deserialize)]
struct BinaryData {
    assembly_result: AssemblyResult,
    graph_nodes: Vec<NodeAttr>,
    graph_edges: Vec<(usize, usize, EdgeAttr)>,
    metadata: Metadata,
}

pub fn write_binary(
    result: &AssemblyResult,
    dag: &AssemblyDag,
    metadata: &Metadata,
    path: &Path,
) -> Result<()> {
    let graph_nodes: Vec<NodeAttr> = dag.node_weights().cloned().collect();
    let graph_edges: Vec<(usize, usize, EdgeAttr)> = dag
        .edge_indices()
        .map(|e| {
            let (s, t) = dag.edge_endpoints(e).unwrap();
            (s.index(), t.index(), dag[e])
        })
        .collect();

    let data = BinaryData {
        assembly_result: result.clone(),
        graph_nodes,
        graph_edges,
        metadata: metadata.clone(),
    };

    let json = serde_json::to_vec(&data)?;
    std::fs::write(path, json)
        .with_context(|| format!("Failed to write binary result: {}", path.display()))?;
    Ok(())
}

pub fn read_binary(path: &Path) -> Result<(AssemblyResult, AssemblyDag, Metadata)> {
    let data = std::fs::read(path)
        .with_context(|| format!("Failed to read binary result: {}", path.display()))?;
    let binary: BinaryData = serde_json::from_slice(&data)?;

    let mut dag = AssemblyDag::new();
    let node_indices: Vec<_> = binary
        .graph_nodes
        .into_iter()
        .map(|attr| dag.add_node(attr))
        .collect();

    for (source, target, attr) in binary.graph_edges {
        dag.add_edge(node_indices[source], node_indices[target], attr);
    }

    Ok((binary.assembly_result, dag, binary.metadata))
}
