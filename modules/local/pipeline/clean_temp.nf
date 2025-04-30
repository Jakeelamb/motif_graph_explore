process clean_temp {
    tag "$db_out"
    
    input:
    path db_out
    
    script:
    """
    # Initialize Conda
    source /nfs/home/jlamb/bin/miniconda3/etc/profile.d/conda.sh
    
    # Activate the motif environment
    conda activate motif
    
    # Clean up temporary files
    rm -rf work/
    rm -f .nextflow.log*
    """
} 