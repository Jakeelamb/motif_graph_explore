process index_genome {
    tag "$genome"
    publishDir "${params.outdir}/genomes", mode: 'copy'
    
    input:
    path genome
    
    output:
    path "${genome}_indexed", emit: indexed_genome
    
    script:
    """
    # No need for manual conda initialization - handled by the config
    
    echo "[INDEX_GENOME] Starting genome indexing for: ${genome} at \$(date)"
    mkdir -p ${genome}_indexed
    
    # Find all genome FASTA files in the ncbi_dataset subdirectory structure
    echo "[INDEX_GENOME] Searching for FASTA files in ${genome}/ncbi_dataset"
    find ${genome}/ncbi_dataset -type f -name "*.fna" > genome_files.txt
    
    # Check if any FASTA files were found
    if [ ! -s genome_files.txt ]; then
        echo "[INDEX_GENOME] ERROR: No FASTA files found in ${genome}/ncbi_dataset"
        exit 1
    else
        echo "[INDEX_GENOME] Found FASTA files: \$(cat genome_files.txt | wc -l) files"
        echo "[INDEX_GENOME] FASTA files: \$(cat genome_files.txt)"
    fi
    
    # Combine all genome files into a single FASTA
    echo "[INDEX_GENOME] Combining FASTA files into combined_genome.fasta"
    cat \$(cat genome_files.txt) > combined_genome.fasta
    
    # Check if the combined file is empty
    if [ ! -s combined_genome.fasta ]; then
        echo "[INDEX_GENOME] ERROR: Combined FASTA file is empty"
        exit 1
    else
        echo "[INDEX_GENOME] Combined FASTA file created successfully, size: \$(ls -lh combined_genome.fasta | awk '{print \\\$5}')"
    fi
    
    # Create BEA-MEM2 index for the genome
    echo "[INDEX_GENOME] Creating BWA-MEM2 index for combined_genome.fasta"
    bwa-mem2 index combined_genome.fasta
    if [ \$? -eq 0 ]; then
        echo "[INDEX_GENOME] BWA-MEM2 index creation successful"
    else
        echo "[INDEX_GENOME] ERROR: BWA-MEM2 index creation failed"
    fi
    
    # Create samtools index
    echo "[INDEX_GENOME] Creating samtools index for combined_genome.fasta"
    samtools faidx combined_genome.fasta
    if [ \$? -eq 0 ]; then
        echo "[INDEX_GENOME] Samtools index creation successful"
    else
        echo "[INDEX_GENOME] ERROR: Samtools index creation failed"
    fi
    
    # Move all files to the indexed directory
    echo "[INDEX_GENOME] Moving indexed files to ${genome}_indexed/"
    mv combined_genome.fasta* ${genome}_indexed/
    
    # Create a completion flag file
    echo "[INDEX_GENOME] Genome indexing completed on \$(date)"
    echo "Genome indexing completed on \$(date)" > ${genome}_indexed/indexing_complete.txt
    """
}

