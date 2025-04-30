process append_to_db {
    tag "${streme_results}_${trf_results}"
    publishDir "${params.outdir}/db", mode: 'copy'
    
    input:
    path streme_results
    path trf_results
    
    output:
    path "motifs.db", emit: motif_db
    
    script:
    """
    # Initialize Conda
    source /home/jake/miniconda3/etc/profile.d/conda.sh
    
    # Activate the motif environment
    conda activate motif
    
    # Create a Python script to update the motif database
    cat > update_motif_db.py << 'EOF'
#!/usr/bin/env python3
import sys
import os
import pandas as pd
import sqlite3
import glob

def create_motif_db(db_path):
    '''Create a new SQLite database with tables for motifs and metadata'''
    conn = sqlite3.connect(db_path)
    
    # Create tables
    conn.execute('''
    CREATE TABLE IF NOT EXISTS motifs (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        motif_id TEXT,
        consensus TEXT,
        width INTEGER,
        source TEXT,
        species TEXT,
        p_value REAL,
        e_value REAL,
        sites INTEGER,
        period INTEGER,
        copies REAL,
        score INTEGER,
        sequence TEXT,
        start INTEGER,
        end INTEGER,
        percent_match REAL
    )
    ''')
    
    conn.execute('''
    CREATE TABLE IF NOT EXISTS metadata (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        species TEXT,
        genome_size INTEGER,
        masked_percent REAL,
        window_size INTEGER,
        window_count INTEGER,
        filtered_window_count INTEGER,
        streme_motif_count INTEGER,
        trf_motif_count INTEGER,
        processing_date TEXT
    )
    ''')
    
    conn.commit()
    conn.close()

def update_motif_db(streme_dir, trf_dir, db_path):
    '''Update the motif database with STREME and TRF results'''
    # Create database if it doesn't exist
    if not os.path.exists(db_path):
        create_motif_db(db_path)
    
    conn = sqlite3.connect(db_path)
    
    # Load STREME motifs
    streme_tsv = glob.glob(f"{streme_dir}/parsed_streme.tsv")
    if streme_tsv and os.path.exists(streme_tsv[0]) and os.path.getsize(streme_tsv[0]) > 0:
        streme_df = pd.read_csv(streme_tsv[0], sep='\\t')
        streme_motif_count = len(streme_df)
        print(f"Loaded {streme_motif_count} motifs from STREME")
        
        # Get species name from directory structure
        species = os.path.basename(streme_dir).split('_')[0]
        streme_df['species'] = species
        
        # Insert STREME motifs into database
        streme_df.to_sql('motifs', conn, if_exists='append', index=False)
    else:
        streme_motif_count = 0
        print("No STREME motifs to load")
    
    # Load TRF motifs
    trf_tsv = glob.glob(f"{trf_dir}/parsed_trf.tsv")
    if trf_tsv and os.path.exists(trf_tsv[0]) and os.path.getsize(trf_tsv[0]) > 0:
        trf_df = pd.read_csv(trf_tsv[0], sep='\\t')
        trf_motif_count = len(trf_df)
        print(f"Loaded {trf_motif_count} motifs from TRF")
        
        # Get species name from directory structure
        species = os.path.basename(trf_dir).split('_')[0]
        trf_df['species'] = species
        
        # Insert TRF motifs into database
        trf_df.to_sql('motifs', conn, if_exists='append', index=False)
    else:
        trf_motif_count = 0
        print("No TRF motifs to load")
    
    # Update metadata
    import datetime
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Read stats files if available
    window_stats_file = glob.glob(f"{streme_dir.replace('_streme', '_windows')}/window_stats.txt")
    window_count = 0
    filtered_window_count = 0
    
    if window_stats_file and os.path.exists(window_stats_file[0]):
        with open(window_stats_file[0], 'r') as f:
            for line in f:
                if "Total windows:" in line:
                    try:
                        window_count = int(line.split(':')[1].strip())
                    except ValueError:
                        pass
                elif "Windows after filtering" in line:
                    try:
                        filtered_window_count = int(line.split(':')[1].strip())
                    except ValueError:
                        pass
    
    # Get species name from one of the result directories
    species = os.path.basename(streme_dir).split('_')[0] or os.path.basename(trf_dir).split('_')[0]
    
    # Check if metadata already exists for this species
    cursor = conn.execute("SELECT id FROM metadata WHERE species = ?", (species,))
    metadata_exists = cursor.fetchone() is not None
    
    if metadata_exists:
        conn.execute('''
        UPDATE metadata SET
        window_count = ?,
        filtered_window_count = ?,
        streme_motif_count = ?,
        trf_motif_count = ?,
        processing_date = ?
        WHERE species = ?
        ''', (window_count, filtered_window_count, streme_motif_count, trf_motif_count, now, species))
    else:
        conn.execute('''
        INSERT INTO metadata
        (species, window_count, filtered_window_count, streme_motif_count, trf_motif_count, processing_date)
        VALUES (?, ?, ?, ?, ?, ?)
        ''', (species, window_count, filtered_window_count, streme_motif_count, trf_motif_count, now))
    
    conn.commit()
    conn.close()
    print(f"Updated motif database with {streme_motif_count} STREME motifs and {trf_motif_count} TRF motifs")

if __name__ == "__main__":
    streme_dir = sys.argv[1]
    trf_dir = sys.argv[2]
    db_path = sys.argv[3]
    
    update_motif_db(streme_dir, trf_dir, db_path)
EOF
    
    # Make the script executable
    chmod +x update_motif_db.py
    
    # Run the database update script
    ./update_motif_db.py ${streme_results} ${trf_results} motifs.db
    """
}