use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;

mod assembly;
mod encode;
mod fasta;
mod graph;
mod kmer;
mod output;

#[derive(Parser)]
#[command(name = "assembly_theory")]
#[command(about = "Compute Assembly Index on genomic DNA sequences")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Run the full assembly index pipeline
    Build {
        /// Input FASTA file
        #[arg(short, long)]
        input: PathBuf,

        /// Maximum k-mer length (default: adaptive)
        #[arg(long)]
        max_k: Option<usize>,

        /// Number of threads (default: all available)
        #[arg(short, long)]
        threads: Option<usize>,

        /// Output directory
        #[arg(short, long, default_value = "output")]
        output: PathBuf,
    },

    /// Export results to various formats
    Export {
        /// Input binary result file
        #[arg(short, long)]
        input: PathBuf,

        /// Output format
        #[arg(short, long, value_enum)]
        format: output::ExportFormat,

        /// Output file path
        #[arg(short, long)]
        output: PathBuf,
    },

    /// Print summary statistics
    Info {
        /// Input binary result file
        #[arg(short, long)]
        input: PathBuf,
    },
}

fn main() -> Result<()> {
    tracing_subscriber::fmt::init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Build {
            input,
            max_k,
            threads,
            output: output_dir,
        } => cmd_build(input, max_k, threads, output_dir),
        Commands::Export {
            input,
            format,
            output: output_path,
        } => cmd_export(input, format, output_path),
        Commands::Info { input } => cmd_info(input),
    }
}

fn cmd_build(
    input: PathBuf,
    max_k: Option<usize>,
    threads: Option<usize>,
    output_dir: PathBuf,
) -> Result<()> {
    use std::time::Instant;

    if let Some(t) = threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .ok();
    }

    let start = Instant::now();

    // 1. Read genome (mmap — no data copied)
    tracing::info!("Reading FASTA: {}", input.display());
    let reader = fasta::FastaReader::from_path(&input)?;
    let genome_length = reader.genome_length();
    tracing::info!(
        "Loaded {} sequences, {} total bases",
        reader.regions().len(),
        genome_length
    );

    // 2. Determine max k
    let max_k = max_k.unwrap_or_else(|| kmer::find_max_k(genome_length));
    tracing::info!("Using max_k = {}", max_k);

    // 3. Count k-mers and build vocabulary (streams from mmap)
    tracing::info!("Counting k-mers...");
    let vocabulary = kmer::build_vocabulary(&reader, max_k);
    tracing::info!("Vocabulary size: {} k-mers", vocabulary.len());

    // 4. Compute assembly index
    tracing::info!("Computing assembly index...");
    let assembly_result = assembly::compute_assembly_index(&vocabulary, max_k);
    tracing::info!(
        "Max assembly index: {}",
        assembly_result.max_assembly_index()
    );

    // 5. Build graph
    tracing::info!("Building assembly DAG...");
    let dag = graph::build_dag(&assembly_result);
    let stats = graph::dag_stats(&dag);
    tracing::info!(
        "DAG: {} nodes, {} edges, max depth {}",
        stats.node_count,
        stats.edge_count,
        stats.max_depth
    );

    // 6. Export all formats
    std::fs::create_dir_all(&output_dir)?;

    let dot_path = output_dir.join("assembly.dot");
    output::write_dot(&dag, &dot_path)?;
    tracing::info!("Wrote DOT: {}", dot_path.display());

    let graphml_path = output_dir.join("assembly.graphml");
    output::write_graphml(&dag, &graphml_path)?;
    tracing::info!("Wrote GraphML: {}", graphml_path.display());

    let csv_path = output_dir.join("metrics.csv");
    output::write_csv(&assembly_result, &csv_path)?;
    tracing::info!("Wrote CSV: {}", csv_path.display());

    let metadata = output::Metadata {
        genome_file: input.display().to_string(),
        genome_length,
        max_k,
        vocab_size: vocabulary.len(),
        max_assembly_index: assembly_result.max_assembly_index(),
        node_count: stats.node_count,
        edge_count: stats.edge_count,
        elapsed_secs: start.elapsed().as_secs_f64(),
    };

    let json_path = output_dir.join("metadata.json");
    output::write_json(&metadata, &json_path)?;
    tracing::info!("Wrote metadata: {}", json_path.display());

    let bin_path = output_dir.join("result.bin");
    output::write_binary(&assembly_result, &dag, &metadata, &bin_path)?;
    tracing::info!("Wrote binary: {}", bin_path.display());

    tracing::info!("Done in {:.2}s", start.elapsed().as_secs_f64());
    tracing::info!(
        "Render: dot -Tsvg {}/assembly.dot -o assembly.svg",
        output_dir.display()
    );
    Ok(())
}

fn cmd_export(input: PathBuf, format: output::ExportFormat, output_path: PathBuf) -> Result<()> {
    let (assembly_result, dag, metadata) = output::read_binary(&input)?;

    match format {
        output::ExportFormat::Graphml => output::write_graphml(&dag, &output_path)?,
        output::ExportFormat::Dot => output::write_dot(&dag, &output_path)?,
        output::ExportFormat::Csv => output::write_csv(&assembly_result, &output_path)?,
        output::ExportFormat::Json => output::write_json(&metadata, &output_path)?,
    }

    tracing::info!("Exported to {}", output_path.display());
    Ok(())
}

fn cmd_info(input: PathBuf) -> Result<()> {
    let (_assembly_result, _dag, metadata) = output::read_binary(&input)?;

    println!("Assembly Theory Results");
    println!("=======================");
    println!("Genome file:          {}", metadata.genome_file);
    println!("Genome length:        {} bp", metadata.genome_length);
    println!("Max k:                {}", metadata.max_k);
    println!("Vocabulary size:      {}", metadata.vocab_size);
    println!("Max assembly index:   {}", metadata.max_assembly_index);
    println!("Graph nodes:          {}", metadata.node_count);
    println!("Graph edges:          {}", metadata.edge_count);
    println!("Computation time:     {:.2}s", metadata.elapsed_secs);

    Ok(())
}
