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
    source /home/jake/miniconda3/etc/profile.d/conda.sh
    
    # Activate the motif environment
    conda activate motif
    
    mkdir -p ${masked_genome}_windows
    
    # Get the masked genome and index files using absolute paths
    MASKED_FASTA="/home/jake/Projects/motif_graph_explore/results/genomes/${masked_genome}/genome.masked.fasta"
    INDEX_FILE="/home/jake/Projects/motif_graph_explore/results/genomes/${masked_genome}/genome.masked.fasta.fai"
    
    # Check if files exist
    if [ ! -f "\${MASKED_FASTA}" ]; then
        echo "Error: Masked FASTA file not found at \${MASKED_FASTA}" >&2
        exit 1
    fi
    if [ ! -f "\${INDEX_FILE}" ]; then
        echo "Error: Index file not found at \${INDEX_FILE}" >&2
        exit 1
    fi
    
    # Generate 10kb windows using bedtools
    bedtools makewindows -g \${INDEX_FILE} -w 10000 > ${masked_genome}_windows.bed
    
    # Extract sequences for each window
    bedtools getfasta -fi \${MASKED_FASTA} -bed ${masked_genome}_windows.bed > ${masked_genome}_windows.fa
    
    # Create a Python script to filter windows with >50% softmasked bases
    cat > filter_softmasked.py << 'EOF'
#!/usr/bin/env python3
import sys
import os

def count_lowercase(seq):
    return sum(1 for c in seq if c.islower())

def main():
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    threshold = 0.5  # 50% threshold
    
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        header = None
        seq = ""
        
        for line in fin:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Process previous sequence if it exists
                if header and seq:
                    lowercase_ratio = count_lowercase(seq) / len(seq) if len(seq) > 0 else 1.0
                    if lowercase_ratio <= threshold:
                        fout.write(header + '\\n')
                        fout.write(seq + '\\n')
                
                # Start new sequence
                header = line
                seq = ""
            else:
                seq += line
        
        # Process the last sequence
        if header and seq:
            lowercase_ratio = count_lowercase(seq) / len(seq) if len(seq) > 0 else 1.0
            if lowercase_ratio <= threshold:
                fout.write(header + '\\n')
                fout.write(seq + '\\n')

if __name__ == "__main__":
    main()
EOF
    
    # Make the script executable
    chmod +x filter_softmasked.py
    
    # Run the filtering script
    ./filter_softmasked.py ${masked_genome}_windows.fa ${masked_genome}_windows/filtered_windows.fa
    
    # Generate statistics
    echo "Window statistics for ${masked_genome}" > ${masked_genome}_windows/window_stats.txt
    echo "Total windows: \$(grep -c '>' ${masked_genome}_windows.fa)" >> ${masked_genome}_windows/window_stats.txt
    echo "Windows after filtering (<=50% softmasked): \$(grep -c '>' ${masked_genome}_windows/filtered_windows.fa)" >> ${masked_genome}_windows/window_stats.txt
    echo "Processing date: \$(date)" >> ${masked_genome}_windows/window_stats.txt
    
    # Copy the BED file to the output directory
    cp ${masked_genome}_windows.bed ${masked_genome}_windows/
    """
}