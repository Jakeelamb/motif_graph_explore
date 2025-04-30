process softmask_genome {
    tag "$indexed_genome"
    publishDir "${params.outdir}/genomes", mode: 'copy'
    
    input:
    path indexed_genome
    
    output:
    path "${indexed_genome}_masked", emit: masked_genome
    
    script:
    """
    # Initialize Conda
    source /nfs/home/jlamb/bin/miniconda3/etc/profile.d/conda.sh
    
    # Activate the motif environment
    conda activate motif
    
    echo "[SOFTMASK_GENOME] Starting genome softmasking for: ${indexed_genome} at \\\$(date)"
    mkdir -p ${indexed_genome}_masked
    
    # Copy the combined genome FASTA file
    echo "[SOFTMASK_GENOME] Copying combined genome FASTA file from ${indexed_genome}/combined_genome.fasta"
    cp ${indexed_genome}/combined_genome.fasta ./genome_to_mask.fasta
    if [ \$? -eq 0 ]; then
        echo "[SOFTMASK_GENOME] FASTA file copied successfully"
    else
        echo "[SOFTMASK_GENOME] ERROR: Failed to copy FASTA file"
    fi
    
    # Skip RepeatMasker for local development to prevent crash
    echo "[SOFTMASK_GENOME] Skipping RepeatMasker softmasking for local run"
    cp ./genome_to_mask.fasta genome_to_mask.fasta.masked
    
    # Move the masked genome and related files to the output directory
    echo "[SOFTMASK_GENOME] Moving masked genome to ${indexed_genome}_masked/genome.masked.fasta"
    mv genome_to_mask.fasta.masked ${indexed_genome}_masked/genome.masked.fasta
    echo "[SOFTMASK_GENOME] Copying index file to ${indexed_genome}_masked/genome.masked.fasta.fai"
    cp ${indexed_genome}/combined_genome.fasta.fai ${indexed_genome}_masked/genome.masked.fasta.fai
    
    # If there's a RepeatMasker summary file, copy it for logging
    if [ -f genome_to_mask.fasta.tbl ]; then
        echo "[SOFTMASK_GENOME] Copying RepeatMasker summary file to ${indexed_genome}_masked/"
        cp genome_to_mask.fasta.tbl ${indexed_genome}_masked/
    else
        echo "[SOFTMASK_GENOME] No RepeatMasker summary file found"
    fi
    
    # Create a completion flag file
    echo "[SOFTMASK_GENOME] RepeatMasker skipped on \\\$(date)"
    echo "RepeatMasker skipped on \\\$(date)" > ${indexed_genome}_masked/masking_complete.txt
    """
}

