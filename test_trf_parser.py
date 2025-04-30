#!/usr/bin/env python3
import os
import sys
import re
import pandas as pd
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='[TRF_TEST] %(levelname)s: %(message)s')

def parse_trf_dat(trf_dat_path, output_tsv):
    '''Parse TRF .dat file and convert to TSV with motif information'''
    
    # Initialize list to store tandem repeat data
    repeats = []
    sequence_name = None # To store the current sequence name

    try:
        with open(trf_dat_path, 'r') as f:
            lines = f.readlines()
        logging.info(f"Read {len(lines)} lines from {trf_dat_path}")
        logging.info(f"First 10 lines of file:")
        for i, line in enumerate(lines[:10]):
            logging.info(f"Line {i+1}: {line.strip()}")

    except FileNotFoundError:
        logging.error(f"Input file not found: {trf_dat_path}")
        return None
    except Exception as e:
        logging.error(f"Error reading file {trf_dat_path}: {e}")
        return None

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

        # Log raw lines for debugging (first 20 data-like lines)
        if len(lines) > 5 and i < 50 and line and not line.startswith("Sequence:") and len(line.split()) > 5:
            fields = line.split()
            found_data_lines += 1
            
            # Check if this looks like a data line
            if len(fields) >= 12 and fields[0].isdigit() and fields[1].isdigit():
                logging.info(f"Data line {found_data_lines} at line {i+1}: '{line}'")
                logging.info(f"Field count: {len(fields)}")
                if found_data_lines <= 3:
                    for j, field in enumerate(fields):
                        logging.info(f"   Field {j}: '{field}'")
                field_indices.append(len(fields))

        # Attempt to parse data lines (assuming space delimited)
        # Format: Start End Period Copies ConsensusSize PercentMatches PercentIndels Score A C G T Entropy ConsensusSequence RepeatSequence
        fields = line.split() # Simple space split

        # Basic check for a plausible data line (needs at least 13 fields typically)
        if len(fields) >= 13 and fields[0].isdigit() and fields[1].isdigit():
            try:
                # TRF-dat format typically has consensus sequence at index 13
                # But we should check field counts to confirm
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
            
            return df
        except Exception as e:
            logging.error(f"Error writing TSV file {output_tsv}: {e}")
            return None
    else:
        logging.warning(f"No tandem repeats found or parsed from {trf_dat_path}")
        return None

if __name__ == "__main__":
    # File to parse
    trf_dat_path = "GCF_001040885.1_S_ratti_ED321_genomic.fna.2.7.7.80.10.50.500.dat"
    output_tsv = "test_parsed_trf.tsv"
    
    logging.info(f"Testing TRF parser with file: {trf_dat_path}")
    
    if not os.path.exists(trf_dat_path):
        logging.error(f"Error: TRF .dat file not found at {trf_dat_path}")
        sys.exit(1)
    elif os.path.getsize(trf_dat_path) == 0:
        logging.warning(f"TRF .dat file is empty: {trf_dat_path}")
        sys.exit(1)
    else:
        logging.info(f"TRF .dat file exists and has size: {os.path.getsize(trf_dat_path)} bytes")
        df = parse_trf_dat(trf_dat_path, output_tsv)
        
        if df is not None and not df.empty:
            # Check if output file was created successfully
            if os.path.exists(output_tsv):
                logging.info(f"Output file created successfully: {output_tsv}")
                logging.info(f"Output file size: {os.path.getsize(output_tsv)} bytes")
                
                # Print first few rows
                logging.info("First 5 rows of parsed data:")
                print(df.head(5).to_string())
            else:
                logging.error(f"Output file was not created: {output_tsv}")
        else:
            logging.error("Parsing failed or no results were found") 