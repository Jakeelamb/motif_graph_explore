#!/bin/bash

# Load necessary modules or activate conda environment
source /nfs/home/jlamb/bin/miniconda3/bin/activate motif

# Run the Nextflow pipeline
nextflow run main.nf -profile slurm -resume

# To run with high memory profile if needed:
# nextflow run main.nf -profile highmem -resume 