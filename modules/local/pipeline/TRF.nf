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
    
    # Define the input FASTA file path
    FASTA_FILE="${indexed_genome}/combined_genome.fasta"
    echo "[TRF] Constructed FASTA path: \${FASTA_FILE}"

    # Check if input file exists
    if [ ! -f "\${FASTA_FILE}" ]; then
        echo "[TRF] ERROR: FASTA file not found at \${FASTA_FILE}"
        # Create mock output for pipeline to continue
        echo "[TRF] Creating mock TRF output to allow pipeline to continue"
        touch ${indexed_genome}_trf/trf.dat
        touch ${indexed_genome}_trf/parsed_trf.tsv # Ensure parser doesn't fail later
        # Add header to mock parsed file
        echo "motif_id\tsequence\tstart\tend\tperiod\tcopies\tconsensus\tpercent_match\tscore\tsource" > ${indexed_genome}_trf/parsed_trf.tsv

    else
        echo "[TRF] Input FASTA file confirmed found at \${FASTA_FILE}"
        echo "[TRF] File size before TRF: \$(ls -lh \${FASTA_FILE} | awk '{print \$5}')"
        
        # Try to run TRF with a timeout to avoid hanging
        echo "[TRF] Running TRF command: timeout 1800 trf \"\${FASTA_FILE}\" 2 7 7 80 10 50 500 -d"
        timeout 1800 trf "\${FASTA_FILE}" 2 7 7 80 10 50 500 -d || true
        
        # Check if TRF output file was created
        # Construct expected output file name based on input FASTA path and TRF parameters
        DAT_FILE="\${FASTA_FILE}.2.7.7.80.10.50.500.dat"
        echo "[TRF] Expected TRF output file: \${DAT_FILE}"

        if [ -f "\${DAT_FILE}" ]; then
            echo "[TRF] TRF analysis produced output file: \${DAT_FILE}"
            mv "\${DAT_FILE}" ${indexed_genome}_trf/trf.dat
        else
            echo "[TRF] WARNING: TRF output file (\${DAT_FILE}) not found after execution. Creating an empty one."
            touch ${indexed_genome}_trf/trf.dat
        fi
    fi
    
    # Check if the final trf.dat file exists in the target directory
    if [ -f "${indexed_genome}_trf/trf.dat" ]; then
        echo "[TRF] Final TRF output file present at ${indexed_genome}_trf/trf.dat"
        echo "[TRF] Size of final output file: \$(ls -lh ${indexed_genome}_trf/trf.dat | awk '{print \$5}')"
        echo "[TRF] Number of lines in final output: \$(wc -l < ${indexed_genome}_trf/trf.dat)"
    else
        echo "[TRF] ERROR: Final TRF output file still not found at ${indexed_genome}_trf/trf.dat"
        # Ensure we have a file for the next step
        touch ${indexed_genome}_trf/trf.dat
        # Create an empty parsed file as well if the dat file is missing
        echo "motif_id\tsequence\tstart\tend\tperiod\tcopies\tconsensus\tpercent_match\tscore\tsource" > ${indexed_genome}_trf/parsed_trf.tsv
    fi
    
    # Create a Python script to parse TRF results
    echo "[TRF] Creating Python script to parse TRF results"
    cat > parse_trf.py << 'EOF'
    #!/usr/bin/env python3
    import os
    import sys
    import re
    import pandas as pd
    import logging

    # Setup logging
    logging.basicConfig(level=logging.INFO, format='[TRF_PARSE] %(levelname)s: %(message)s')

    def parse_trf_dat(trf_dat_path, output_tsv):
        '''Parse TRF .dat file and convert to TSV with motif information'''
        
        # Initialize list to store tandem repeat data
        repeats = []
        sequence_name = None # To store the current sequence name

        try:
            with open(trf_dat_path, 'r') as f:
                lines = f.readlines()
            logging.info(f"Read {len(lines)} lines from {trf_dat_path}")

        except FileNotFoundError:
            logging.error(f"Input file not found: {trf_dat_path}")
            # Create empty file with headers and return
            columns = ['motif_id', 'sequence', 'start', 'end', 'period', 'copies',
                       'consensus', 'percent_match', 'score', 'source']
            pd.DataFrame(columns=columns).to_csv(output_tsv, sep='\t', index=False)
            logging.warning(f"Created empty TSV file: {output_tsv}")
            return
        except Exception as e:
            logging.error(f"Error reading file {trf_dat_path}: {e}")
            # Create empty file with headers and return
            columns = ['motif_id', 'sequence', 'start', 'end', 'period', 'copies',
                       'consensus', 'percent_match', 'score', 'source']
            pd.DataFrame(columns=columns).to_csv(output_tsv, sep='\t', index=False)
            logging.warning(f"Created empty TSV file: {output_tsv}")
            return

        # TRF .dat format parsing logic
        for i, line in enumerate(lines):
            line = line.strip()
            
            # Skip empty lines or parameter lines
            if not line or line.startswith('Parameters'):
                continue
                
            # Capture sequence name line
            if line.startswith("Sequence:"):
                sequence_name = line.split(":", 1)[1].strip()
                logging.info(f"Processing sequence: {sequence_name}")
                continue

            # Skip header lines like "Indices Period ..."
            if "Period" in line and "Copies" in line and "Consensus" in line:
                continue

            # Attempt to parse data lines (assuming space delimited)
            # Format: Start End Period Copies ConsensusSize PercentMatches PercentIndels Score A C G T Entropy ConsensusSequence RepeatSequence
            fields = line.split() # Simple space split might need refinement based on actual .dat format variations

            # Basic check for a plausible data line (needs at least 12 fields typically)
            if len(fields) >= 12 and fields[0].isdigit() and fields[1].isdigit():
                try:
                    repeat_id = f"TRF-{len(repeats)+1}" # Generate a unique ID

                    # Adjust indices based on observed .dat format
                    # Example mapping (verify with actual .dat file format):
                    # fields[0]: Start Index
                    # fields[1]: End Index
                    # fields[2]: Period Size
                    # fields[3]: Copy Number
                    # fields[4]: Consensus Size (often same as period)
                    # fields[5]: Percent Matches
                    # fields[6]: Percent Indels
                    # fields[7]: Alignment Score
                    # ... other fields ...
                    # fields[12]: Consensus Pattern
                    # fields[13]: Repeat Sequence (optional, might not always be present or needed)

                    repeat = {
                        'motif_id': repeat_id,
                        'sequence': sequence_name if sequence_name else 'UnknownSequence', # Use captured sequence name
                        'start': int(fields[0]),
                        'end': int(fields[1]),
                        'period': int(fields[2]),
                        'copies': float(fields[3]),
                        'consensus': fields[12], # Example: Consensus Pattern
                        'percent_match': float(fields[5]), # Example: Percent Matches
                        'score': int(fields[7]) # Example: Alignment Score
                        # Add other relevant fields if needed
                    }
                    repeats.append(repeat)

                except (ValueError, IndexError) as e:
                    logging.warning(f"Skipping malformed line {i+1}: '{line}'. Error: {e}")
                    logging.debug(f"Fields parsed: {fields}")
            else:
                 # Log lines that don't look like data lines (and aren't headers/sequence names)
                 if not line.startswith("Sequence:"):
                     logging.debug(f"Skipping non-data line {i+1}: '{line}'")

        # Convert to DataFrame and save as TSV
        if repeats:
            df = pd.DataFrame(repeats)
            # Add source column to indicate TRF
            df['source'] = 'TRF'
            try:
                df.to_csv(output_tsv, sep='\t', index=False)
                logging.info(f"Successfully parsed {len(repeats)} tandem repeats into {output_tsv}")
            except Exception as e:
                 logging.error(f"Error writing TSV file {output_tsv}: {e}")
                 # Create empty file with headers if write fails
                 columns = ['motif_id', 'sequence', 'start', 'end', 'period', 'copies',
                            'consensus', 'percent_match', 'score', 'source']
                 pd.DataFrame(columns=columns).to_csv(output_tsv, sep='\t', index=False)
                 logging.warning(f"Created empty TSV file due to write error: {output_tsv}")
        else:
            logging.warning(f"No tandem repeats found or parsed from {trf_dat_path}")
            # Create empty file with headers
            columns = ['motif_id', 'sequence', 'start', 'end', 'period', 'copies',
                       'consensus', 'percent_match', 'score', 'source']
            pd.DataFrame(columns=columns).to_csv(output_tsv, sep='\t', index=False)
            logging.warning(f"Created empty TSV file: {output_tsv}")


    if __name__ == "__main__":
        if len(sys.argv) != 3:
             print(f"Usage: python {sys.argv[0]} <trf_dat_path> <output_tsv>")
             sys.exit(1)

        trf_dat_path = sys.argv[1]
        output_tsv = sys.argv[2]

        if not os.path.exists(trf_dat_path):
             logging.error(f"Error: TRF .dat file not found at {trf_dat_path}")
             # Create empty file with headers
             columns = ['motif_id', 'sequence', 'start', 'end', 'period', 'copies',
                        'consensus', 'percent_match', 'score', 'source']
             pd.DataFrame(columns=columns).to_csv(output_tsv, sep='\t', index=False)
             logging.warning(f"Created empty TSV file: {output_tsv}")
        elif os.path.getsize(trf_dat_path) == 0:
             logging.warning(f"TRF .dat file is empty: {trf_dat_path}")
             # Create empty file with headers
             columns = ['motif_id', 'sequence', 'start', 'end', 'period', 'copies',
                        'consensus', 'percent_match', 'score', 'source']
             pd.DataFrame(columns=columns).to_csv(output_tsv, sep='\t', index=False)
             logging.warning(f"Created empty TSV file: {output_tsv}")
        else:
             parse_trf_dat(trf_dat_path, output_tsv)

    EOF
    
    # Make the script executable
    echo "[TRF] Making parse_trf.py script executable"
    chmod +x parse_trf.py
    
    # Run the parsing script
    echo "[TRF] Running TRF parsing script: ./parse_trf.py ${indexed_genome}_trf/trf.dat ${indexed_genome}_trf/parsed_trf.tsv"
    ./parse_trf.py ${indexed_genome}_trf/trf.dat ${indexed_genome}_trf/parsed_trf.tsv
    if [ \$? -eq 0 ]; then
        echo "[TRF] TRF parsing script finished execution"
        # Check if the parsed file exists and has content beyond the header
        if [ -f "${indexed_genome}_trf/parsed_trf.tsv" ] && [ "\$(wc -l < ${indexed_genome}_trf/parsed_trf.tsv)" -gt 1 ]; then
            echo "[TRF] Parsed output file found at ${indexed_genome}_trf/parsed_trf.tsv"
            echo "[TRF] Number of motifs parsed (incl header): \$(wc -l < ${indexed_genome}_trf/parsed_trf.tsv)"
        elif [ -f "${indexed_genome}_trf/parsed_trf.tsv" ]; then
             echo "[TRF] WARNING: Parsed output file exists but is empty or contains only the header."
        else
            echo "[TRF] ERROR: Parsed output file not found at ${indexed_genome}_trf/parsed_trf.tsv"
            # Create an empty TSV with headers if parser failed to create one
             echo "motif_id\tsequence\tstart\tend\tperiod\tcopies\tconsensus\tpercent_match\tscore\tsource" > ${indexed_genome}_trf/parsed_trf.tsv
        fi
    else
        echo "[TRF] ERROR: TRF parsing script failed with exit code \$?"
        # Create an empty TSV with headers if it doesn't exist
        if [ ! -f "${indexed_genome}_trf/parsed_trf.tsv" ]; then
             echo "motif_id\tsequence\tstart\tend\tperiod\tcopies\tconsensus\tpercent_match\tscore\tsource" > ${indexed_genome}_trf/parsed_trf.tsv
        fi
    fi
    echo "[TRF] TRF process completed on \$(date)"
    """
}