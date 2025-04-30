process index_genome {
    tag "$genome"
    publishDir "${params.outdir}/genomes", mode: 'copy'
    
    input:
    path genome
    
    output:
    path "${genome}_indexed", emit: indexed_genome
    
    script:
    """
    # Initialize Conda
    source /home/jake/miniconda3/etc/profile.d/conda.sh
    
    # Activate the motif environment
    conda activate motif
    
    mkdir -p ${genome}_indexed
    
    # Find all genome FASTA files in the ncbi_dataset subdirectory structure
    find ${genome}/ncbi_dataset -type f -name "*.fna" > genome_files.txt
    
    # Check if any FASTA files were found
    if [ ! -s genome_files.txt ]; then
        echo "Error: No FASTA files found in ${genome}/ncbi_dataset" >&2
        exit 1
    fi
    
    # Combine all genome files into a single FASTA
    cat \$(cat genome_files.txt) > combined_genome.fasta
    
    # Check if the combined file is empty
    if [ ! -s combined_genome.fasta ]; then
        echo "Error: Combined FASTA file is empty" >&2
        exit 1
    fi
    
    # Create BEA-MEM2 index for the genome
    bwa-mem2 index combined_genome.fasta
    
    # Create samtools index
    samtools faidx combined_genome.fasta
    
    # Move all files to the indexed directory
    mv combined_genome.fasta* ${genome}_indexed/
    
    # Create a completion flag file
    echo "Genome indexing completed on \$(date)" > ${genome}_indexed/indexing_complete.txt
    """
}

