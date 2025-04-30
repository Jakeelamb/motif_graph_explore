#!/bin/bash -ue
# Initialize Conda
source /home/jake/miniconda3/etc/profile.d/conda.sh

# Activate the motif environment
conda activate motif

echo "[FETCH_GENOME] Starting genome download for accession: GCA_000002985.3 at $(date)"
mkdir -p GCA_000002985.3_genome

# Download genome using NCBI datasets
echo "[FETCH_GENOME] Running datasets download command for GCA_000002985.3"
datasets download genome accession GCA_000002985.3 --filename GCA_000002985.3.zip
if [ $? -eq 0 ]; then
    echo "[FETCH_GENOME] Download successful for GCA_000002985.3"
else
    echo "[FETCH_GENOME] ERROR: Download failed for GCA_000002985.3"
fi

# Unzip the downloaded genome
echo "[FETCH_GENOME] Unzipping genome data for GCA_000002985.3"
unzip GCA_000002985.3.zip -d GCA_000002985.3_genome/
if [ $? -eq 0 ]; then
    echo "[FETCH_GENOME] Unzip successful for GCA_000002985.3"
else
    echo "[FETCH_GENOME] ERROR: Unzip failed for GCA_000002985.3"
fi

# Create a completion flag file
echo "[FETCH_GENOME] Genome download completed on $(date)"
echo "Genome download completed on $(date)" > GCA_000002985.3_genome/download_complete.txt
