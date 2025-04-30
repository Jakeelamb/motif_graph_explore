# Script to import parsed STREME and TRF motif data into motifs.db

import sqlite3
import pandas as pd
import os
import sys

def import_tsv_to_db(tsv_path, source, db_path='motifs.db'):
    if not os.path.exists(tsv_path):
        print(f"File {tsv_path} not found.")
        return 0
    
    df = pd.read_csv(tsv_path, sep='\t')
    if df.empty:
        print(f"No data found in {tsv_path}.")
        return 0
    
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    
    count = 0
    for index, row in df.iterrows():
        cur.execute('''
            INSERT INTO motifs (
                motif_id, consensus, width, source, species, 
                p_value, e_value, sites, period, copies, score, 
                sequence, start, end, percent_match
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            row.get('motif_id', ''),
            row.get('consensus', ''),
            row.get('width', 0),
            source,
            row.get('species', ''),
            row.get('p_value', 0.0),
            row.get('e_value', 0.0),
            row.get('sites', 0),
            row.get('period', 0),
            row.get('copies', 0.0),
            row.get('score', 0),
            row.get('sequence', ''),
            row.get('start', 0),
            row.get('end', 0),
            row.get('percent_match', 0.0)
        ))
        count += 1
    
    conn.commit()
    conn.close()
    print(f"Imported {count} motifs from {tsv_path} into {db_path}.")
    return count

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python import_motifs_to_db.py <path_to_parsed_tsv> [db_path]")
        sys.exit(1)
    
    tsv_path = sys.argv[1]
    db_path = sys.argv[2] if len(sys.argv) > 2 else 'motifs.db'
    source = 'STREME' if 'streme' in tsv_path.lower() else 'TRF'
    
    import_tsv_to_db(tsv_path, source, db_path) 