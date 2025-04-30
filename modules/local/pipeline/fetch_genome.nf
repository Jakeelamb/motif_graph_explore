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
    
    mkdir -p ${accession}_genome
    
    # Download genome using NCBI datasets
    datasets download genome accession ${accession} --filename ${accession}.zip
    
    # Unzip the downloaded genome
    unzip ${accession}.zip -d ${accession}_genome/
    
    # Create a completion flag file
    echo "Genome download completed on \$(date)" > ${accession}_genome/download_complete.txt
    """
} 