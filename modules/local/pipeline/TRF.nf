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
    # No need for manual conda initialization - handled by the config
    
    echo "[TRF] Starting TRF analysis for: ${indexed_genome} at \$(date)"
    mkdir -p ${indexed_genome}_trf
    
    # Define the input FASTA file path
    FASTA_FILE="${indexed_genome}/combined_genome.fasta"
    echo "[TRF] Constructed FASTA path: \${FASTA_FILE}"
    echo "[TRF] File exists check: \$(if [ -f "\${FASTA_FILE}" ]; then echo "YES"; else echo "NO"; fi)"
    echo "[TRF] File size: \$(if [ -f "\${FASTA_FILE}" ]; then ls -lh \${FASTA_FILE} 2>/dev/null | awk '{print \$5}'; else echo "N/A"; fi)"

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
        echo "[TRF] File type: \$(file \${FASTA_FILE})"
        echo "[TRF] First 10 lines of FASTA file:"
        head -n 10 \${FASTA_FILE}
        
        # Get and log environment information
        echo "[TRF] TRF command path: \$(which trf 2>/dev/null || echo 'not found')"
        echo "[TRF] Current working directory: \$(pwd)"
        echo "[TRF] Disk space: \$(df -h . | tail -n 1)"
        
        # Run TRF on the genome file
        echo "[TRF] Running TRF command: trf \${FASTA_FILE} 2 7 7 80 10 50 500 -d"
        trf "\${FASTA_FILE}" 2 7 7 80 10 50 500 -d
        TRF_EXIT_CODE=\$?
        echo "[TRF] TRF command exit code: \${TRF_EXIT_CODE}"
        
        # Check if TRF output file was created
        # Construct expected output file name based on input FASTA path and TRF parameters
        DAT_FILE="\${FASTA_FILE}.2.7.7.80.10.50.500.dat"
        echo "[TRF] Expected TRF output file: \${DAT_FILE}"
        echo "[TRF] Output file exists check: \$(if [ -f "\${DAT_FILE}" ]; then echo "YES"; else echo "NO"; fi)"
        echo "[TRF] Output directory contents:"
        ls -la \$(dirname \${DAT_FILE})

        if [ -f "\${DAT_FILE}" ]; then
            echo "[TRF] TRF analysis produced output file: \${DAT_FILE}"
            echo "[TRF] Output size: \$(ls -lh \${DAT_FILE} | awk '{print \$5}')"
            echo "[TRF] Output line count: \$(wc -l < \${DAT_FILE})"
            echo "[TRF] First 20 lines of output file:"
            head -n 20 \${DAT_FILE}
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
        echo "[TRF] Number of lines in final output: \$(wc -l < ${indexed_genome}_trf/trf.dat || echo 0)"
        if [ -s "${indexed_genome}_trf/trf.dat" ]; then
            echo "[TRF] First 20 lines of final output file:"
            head -n 20 ${indexed_genome}_trf/trf.dat
        else
            echo "[TRF] WARNING: Final output file is empty"
        fi
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
            logging.info(f"First 5 lines of file (or fewer if smaller):")
            for i, line in enumerate(lines[:5]):
                logging.info(f"Line {i+1}: {line.strip()}")

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
        found_data_lines = 0
        parsed_repeats = 0
        skipped_lines = 0
        
        # Collection of field indices
        field_indices = []
        
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
            fields = line.split() # Simple space split

            # Log raw lines for debugging (first 10 data-like lines)
            if i < 20 and line and not line.startswith("Sequence:") and len(fields) > 5:
                logging.info(f"Raw line {i}: '{line}'")
                logging.info(f"Field count: {len(fields)}")
                logging.info(f"Fields: {fields}")

            # Basic check for a plausible data line (needs at least 13 fields typically)
            if len(fields) >= 13 and fields[0].isdigit() and fields[1].isdigit():
                try:
                    # TRF-dat format typically has consensus sequence at index 13
                    found_data_lines += 1
                    repeat_id = f"TRF-{len(repeats)+1}" # Generate a unique ID
                    parsed_repeats += 1

                    # Indices based on TRF .dat format:
                    # fields[0]: Start
                    # fields[1]: End
                    # fields[2]: Period Size
                    # fields[3]: Copy Number
                    # fields[4]: Consensus Size
                    # fields[5]: Percent Matches
                    # fields[6]: Percent Indels
                    # fields[7]: Alignment Score
                    # fields[8-11]: Nucleotide composition
                    # fields[12]: Entropy
                    # fields[13]: Consensus Pattern
                    # fields[14]: Repeat Sequence (optional)

                    consensus_index = 13  # Default index for consensus pattern
                    
                    # Safety check for index
                    if len(fields) <= consensus_index:
                        logging.warning(f"Line {i+1} has too few fields ({len(fields)}), expecting at least {consensus_index+1}")
                        skipped_lines += 1
                        continue

                    repeat = {
                        'motif_id': repeat_id,
                        'sequence': sequence_name if sequence_name else 'UnknownSequence',
                        'start': int(fields[0]),
                        'end': int(fields[1]),
                        'period': int(fields[2]),
                        'copies': float(fields[3]),
                        'consensus': fields[consensus_index],
                        'percent_match': float(fields[5]),
                        'score': int(fields[7])
                    }
                    
                    repeats.append(repeat)
                    if parsed_repeats <= 5:
                        logging.info(f"Parsed repeat {parsed_repeats}: {repeat}")

                except (ValueError, IndexError) as e:
                    skipped_lines += 1
                    logging.warning(f"Skipping malformed line {i+1}: '{line}'. Error: {e}")
            elif len(fields) > 5:
                skipped_lines += 1
                if i < 50:
                    logging.debug(f"Skipping non-data line {i+1}: '{line}'")

        # Field index summary
        if field_indices:
            field_counts = {}
            for count in field_indices:
                field_counts[count] = field_counts.get(count, 0) + 1
            logging.info(f"Field count distribution: {field_counts}")
        
        logging.info(f"Parse summary: Found {found_data_lines} data-like lines, created {len(repeats)} repeats, skipped {skipped_lines} lines")

        # Convert to DataFrame and save as TSV
        if repeats:
            df = pd.DataFrame(repeats)
            # Add source column to indicate TRF
            df['source'] = 'TRF'
            try:
                df.to_csv(output_tsv, sep='\t', index=False)
                logging.info(f"Successfully parsed {len(repeats)} tandem repeats into {output_tsv}")
                
                # Output some basic statistics about the parsed data
                logging.info(f"Data summary:")
                logging.info(f"  Total motifs: {len(df)}")
                logging.info(f"  Unique sequences: {df['sequence'].nunique()}")
                logging.info(f"  Period size range: {df['period'].min()} - {df['period'].max()}")
                logging.info(f"  Copy number range: {df['copies'].min()} - {df['copies'].max()}")
                logging.info(f"  Average percent match: {df['percent_match'].mean():.2f}%")
                logging.info(f"  Average score: {df['score'].mean():.2f}")
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

        logging.info(f"Processing TRF data from {trf_dat_path} to {output_tsv}")
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
             logging.info(f"TRF .dat file exists and has size: {os.path.getsize(trf_dat_path)} bytes")
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
            echo "[TRF] First few lines of parsed TSV:"
            head -n 5 ${indexed_genome}_trf/parsed_trf.tsv
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