// run TRF on the orignal genome fasta file

process TRF {
    tag "$indexed_genome"
    publishDir "${params.outdir}/trf", mode: 'copy'
    
    input:
    path indexed_genome
    
    output:
    path "${indexed_genome}_trf", emit: trf_results
    
    script:
    """
    # Initialize Conda
    source /home/jake/miniconda3/etc/profile.d/conda.sh
    
    # Activate the motif environment
    conda activate motif
    
    mkdir -p ${indexed_genome}_trf
    
    # Run trf on the full genome with a single shell command block using a here-document
    bash << 'EOF'
    FASTA_FILE="/home/jake/Projects/motif_graph_explore/results/genomes/${indexed_genome}/combined_genome.fasta"
    if [ ! -f "\$FASTA_FILE" ]; then
      echo "Error: FASTA file not found at \$FASTA_FILE" >&2
      exit 1
    fi
    trf "\$FASTA_FILE" 2 7 7 80 10 50 500 -d
    mv "\$FASTA_FILE.2.7.7.80.10.50.500.dat" ${indexed_genome}_trf/trf.dat
    EOF
    
    # Create a Python script to parse TRF results
    cat > parse_trf.py << 'EOF'
#!/usr/bin/env python3
import os
import sys
import re
import pandas as pd

def parse_trf_dat(trf_dat_path, output_tsv):
    '''Parse TRF .dat file and convert to TSV with motif information'''
    
    # Initialize list to store tandem repeat data
    repeats = []
    
    with open(trf_dat_path, 'r') as f:
        lines = f.readlines()
    
    # TRF .dat format has a header line followed by data lines
    # Typical data line format: 
    # Sequence  Start  End  Period  Copy#  Consensus  Percent Indels  Score
    
    for line in lines:
        line = line.strip()
        
        # Skip header lines and empty lines
        if not line or line.startswith('Sequence') or line.startswith('Parameters'):
            continue
        
        # Data line
        fields = re.split(r'\\s+', line, maxsplit=8)
        if len(fields) >= 8:
            try:
                repeat_id = f"TRF-{len(repeats)+1}"
                
                repeat = {
                    'motif_id': repeat_id,
                    'sequence': fields[0],
                    'start': int(fields[1]),
                    'end': int(fields[2]),
                    'period': int(fields[3]),
                    'copies': float(fields[4]),
                    'consensus': fields[5],
                    'percent_match': float(re.sub(r'[^\\d.]', '', fields[6])) if fields[6] != '-' else None,
                    'score': int(fields[7])
                }
                
                repeats.append(repeat)
            except (ValueError, IndexError) as e:
                print(f"Error parsing line: {line}")
                print(f"Exception: {e}")
    
    # Convert to DataFrame and save as TSV
    if repeats:
        df = pd.DataFrame(repeats)
        # Add source column to indicate TRF
        df['source'] = 'TRF'
        df.to_csv(output_tsv, sep='\\t', index=False)
        print(f"Parsed {len(repeats)} tandem repeats from TRF output")
    else:
        print("No tandem repeats found in TRF output")
        # Create empty file with headers
        columns = ['motif_id', 'sequence', 'start', 'end', 'period', 'copies', 
                   'consensus', 'percent_match', 'score', 'source']
        pd.DataFrame(columns=columns).to_csv(output_tsv, sep='\\t', index=False)

if __name__ == "__main__":
    trf_dat_path = sys.argv[1]
    output_tsv = sys.argv[2]
    
    if os.path.exists(trf_dat_path):
        parse_trf_dat(trf_dat_path, output_tsv)
    else:
        print(f"Error: TRF .dat file not found at {trf_dat_path}")
        # Create empty file with headers
        columns = ['motif_id', 'sequence', 'start', 'end', 'period', 'copies', 
                   'consensus', 'percent_match', 'score', 'source']
        pd.DataFrame(columns=columns).to_csv(output_tsv, sep='\\t', index=False)
EOF
    
    # Make the script executable
    chmod +x parse_trf.py
    
    # Run the parsing script
    ./parse_trf.py ${indexed_genome}_trf/trf.dat ${indexed_genome}_trf/parsed_trf.tsv
    """
}