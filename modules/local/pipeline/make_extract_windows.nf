process make_extract_windows {
    tag "$masked_genome"
    publishDir "${params.outdir}/windows", mode: 'copy'
    
    input:
    path masked_genome
    
    output:
    path "${masked_genome}_windows", emit: windows
    
    script:
    """
    # Initialize Conda
    source /nfs/home/jlamb/bin/miniconda3/etc/profile.d/conda.sh
    
    # Activate the motif environment
    conda activate motif
    
    echo "[MAKE_EXTRACT_WINDOWS] Starting window creation for: ${masked_genome} at \$(date)"
    mkdir -p ${masked_genome}_windows
    
    # Get the masked genome and index files using absolute paths
    MASKED_FASTA="\$PWD/${params.outdir}/genomes/${masked_genome}/combined_genome.fasta"
    INDEX_FILE="\$PWD/${params.outdir}/genomes/${masked_genome}/combined_genome.fasta.fai"
    
    # Check if files exist
    if [ ! -f "\${MASKED_FASTA}" ]; then
        echo "[MAKE_EXTRACT_WINDOWS] ERROR: Masked FASTA file not found at \${MASKED_FASTA}"
        exit 1
    else
        echo "[MAKE_EXTRACT_WINDOWS] Masked FASTA file found at \${MASKED_FASTA}"
    fi
    if [ ! -f "\${INDEX_FILE}" ]; then
        echo "[MAKE_EXTRACT_WINDOWS] ERROR: Index file not found at \${INDEX_FILE}"
        exit 1
    else
        echo "[MAKE_EXTRACT_WINDOWS] Index file found at \${INDEX_FILE}"
    fi
    
    # Generate 10kb windows using bedtools
    echo "[MAKE_EXTRACT_WINDOWS] Generating 10kb windows using bedtools makewindows"
    bedtools makewindows -g \${INDEX_FILE} -w 10000 > ${masked_genome}_windows.bed
    if [ \$? -eq 0 ]; then
        echo "[MAKE_EXTRACT_WINDOWS] Window generation successful, created ${masked_genome}_windows.bed"
        echo "[MAKE_EXTRACT_WINDOWS] Number of windows generated: \$(wc -l < ${masked_genome}_windows.bed)"
    else
        echo "[MAKE_EXTRACT_WINDOWS] ERROR: Window generation failed"
    fi
    
    # Extract sequences for each window
    echo "[MAKE_EXTRACT_WINDOWS] Extracting sequences for windows using bedtools getfasta"
    bedtools getfasta -fi \${MASKED_FASTA} -bed ${masked_genome}_windows.bed > ${masked_genome}_windows.fa
    if [ \$? -eq 0 ]; then
        echo "[MAKE_EXTRACT_WINDOWS] Sequence extraction successful, created ${masked_genome}_windows.fa"
        echo "[MAKE_EXTRACT_WINDOWS] Number of sequences extracted: \$(grep -c '>' ${masked_genome}_windows.fa)"
    else
        echo "[MAKE_EXTRACT_WINDOWS] ERROR: Sequence extraction failed"
    fi
    
    # Move the unfiltered windows to the output directory
    echo "[MAKE_EXTRACT_WINDOWS] Copying unfiltered windows to output directory"
    mkdir -p ${masked_genome}_windows
    cp ${masked_genome}_windows.fa ${masked_genome}_windows/filtered_windows.fa
    cp ${masked_genome}_windows.bed ${masked_genome}_windows/
    
    # Generate statistics
    echo "[MAKE_EXTRACT_WINDOWS] Generating window statistics"
    echo "Window statistics for ${masked_genome}" > ${masked_genome}_windows/window_stats.txt
    echo "Total windows: \$(grep -c '>' ${masked_genome}_windows.fa)" >> ${masked_genome}_windows/window_stats.txt
    echo "Windows after filtering (<=50% softmasked): \$(grep -c '>' ${masked_genome}_windows/filtered_windows.fa)" >> ${masked_genome}_windows/window_stats.txt
    echo "Processing date: \$(date)" >> ${masked_genome}_windows/window_stats.txt
    echo "[MAKE_EXTRACT_WINDOWS] Window processing completed on \$(date)"
    """
}