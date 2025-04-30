process fetch_genome {
    tag "$accession"
    publishDir "${params.outdir}/genomes", mode: 'copy'
    
    input:
    val accession
    
    output:
    path "${accession}_genome", emit: genome
    
    script:
    """
    # Initialize Conda
    source /home/jake/miniconda3/etc/profile.d/conda.sh
    
    # Activate the motif environment
    conda activate motif
    
    echo "[FETCH_GENOME] Starting genome download for accession: ${accession} at \$(date)"
    mkdir -p ${accession}_genome
    
    # Download genome using NCBI datasets
    echo "[FETCH_GENOME] Running datasets download command for ${accession}"
    datasets download genome accession ${accession} --filename ${accession}.zip
    if [ \$? -eq 0 ]; then
        echo "[FETCH_GENOME] Download successful for ${accession}"
    else
        echo "[FETCH_GENOME] ERROR: Download failed for ${accession}"
    fi
    
    # Unzip the downloaded genome
    echo "[FETCH_GENOME] Unzipping genome data for ${accession}"
    unzip ${accession}.zip -d ${accession}_genome/
    if [ \$? -eq 0 ]; then
        echo "[FETCH_GENOME] Unzip successful for ${accession}"
    else
        echo "[FETCH_GENOME] ERROR: Unzip failed for ${accession}"
    fi
    
    # Create a completion flag file
    echo "[FETCH_GENOME] Genome download completed on \$(date)"
    echo "Genome download completed on \$(date)" > ${accession}_genome/download_complete.txt
    """
} 