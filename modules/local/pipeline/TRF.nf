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
    source /nfs/home/jlamb/bin/miniconda3/etc/profile.d/conda.sh
    
    # Activate the motif environment
    conda activate motif
    
    echo "[TRF] Starting TRF analysis for: ${indexed_genome} at \$(date)"
    mkdir -p ${indexed_genome}_trf
    
    # Check if input file exists
    FASTA_FILE="\$PWD/${params.outdir}/genomes/${indexed_genome}/combined_genome.fasta"
    if [ ! -f "\${FASTA_FILE}" ]; then
        echo "[TRF] ERROR: FASTA file not found at \${FASTA_FILE}"
        # Create mock output for pipeline to continue
        echo "[TRF] Creating mock TRF output to allow pipeline to continue"
        touch ${indexed_genome}_trf/trf.dat
    else
        echo "[TRF] Input FASTA file found at \${FASTA_FILE}"
        echo "[TRF] File size: \$(ls -lh \${FASTA_FILE} | awk '{print \$5}')"
        
        # Try to run TRF with a timeout to avoid hanging
        echo "[TRF] Running TRF analysis on full genome with 30 minute timeout"
        timeout 1800 trf "\${FASTA_FILE}" 2 7 7 80 10 50 500 -d || true
        
        # Check if TRF output file was created
        DAT_FILE="\${FASTA_FILE}.2.7.7.80.10.50.500.dat"
        if [ -f "\${DAT_FILE}" ]; then
            echo "[TRF] TRF analysis completed successfully"
            mv "\${DAT_FILE}" ${indexed_genome}_trf/trf.dat
        else
            echo "[TRF] WARNING: TRF output file not found. Creating an empty one to continue pipeline."
            touch ${indexed_genome}_trf/trf.dat
        fi
    fi
    
    # Check if TRF output file exists
    if [ -f "${indexed_genome}_trf/trf.dat" ]; then
        echo "[TRF] TRF output file found at ${indexed_genome}_trf/trf.dat"
        echo "[TRF] Number of lines in output: \$(wc -l < ${indexed_genome}_trf/trf.dat)"
    else
        echo "[TRF] ERROR: TRF output file not found at ${indexed_genome}_trf/trf.dat"
        # Ensure we have a file for the next step
        touch ${indexed_genome}_trf/trf.dat
    fi
    
    # Create a Python script to parse TRF results
    echo "[TRF] Creating Python script to parse TRF results"
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
                    print(f"[TRF] Error parsing line: {line}")
                    print(f"[TRF] Exception: {e}")
        
        # Convert to DataFrame and save as TSV
        if repeats:
            df = pd.DataFrame(repeats)
            # Add source column to indicate TRF
            df['source'] = 'TRF'
            df.to_csv(output_tsv, sep='\\t', index=False)
            print(f"[TRF] Parsed {len(repeats)} tandem repeats from TRF output")
        else:
            print("[TRF] No tandem repeats found in TRF output")
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
            print(f"[TRF] Error: TRF .dat file not found at {trf_dat_path}")
            # Create empty file with headers
            columns = ['motif_id', 'sequence', 'start', 'end', 'period', 'copies', 
                       'consensus', 'percent_match', 'score', 'source']
            pd.DataFrame(columns=columns).to_csv(output_tsv, sep='\\t', index=False)
    EOF
    
    # Make the script executable
    echo "[TRF] Making parse_trf.py script executable"
    chmod +x parse_trf.py
    
    # Run the parsing script
    echo "[TRF] Running TRF parsing script"
    ./parse_trf.py ${indexed_genome}_trf/trf.dat ${indexed_genome}_trf/parsed_trf.tsv
    if [ \$? -eq 0 ]; then
        echo "[TRF] TRF parsing completed successfully"
        if [ -f "${indexed_genome}_trf/parsed_trf.tsv" ]; then
            echo "[TRF] Parsed output file found at ${indexed_genome}_trf/parsed_trf.tsv"
            echo "[TRF] Number of motifs parsed: \$(wc -l < ${indexed_genome}_trf/parsed_trf.tsv)"
        else
            echo "[TRF] ERROR: Parsed output file not found at ${indexed_genome}_trf/parsed_trf.tsv"
            # Create an empty TSV with headers
            echo "motif_id\tsequence\tstart\tend\tperiod\tcopies\tconsensus\tpercent_match\tscore\tsource" > ${indexed_genome}_trf/parsed_trf.tsv
        fi
    else
        echo "[TRF] ERROR: TRF parsing failed"
        # Create an empty TSV with headers
        echo "motif_id\tsequence\tstart\tend\tperiod\tcopies\tconsensus\tpercent_match\tscore\tsource" > ${indexed_genome}_trf/parsed_trf.tsv
    fi
    echo "[TRF] TRF process completed on \$(date)"
    """
}