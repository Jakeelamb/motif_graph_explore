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
    source /nfs/home/jlamb/bin/miniconda3/etc/profile.d/conda.sh
    
    # Activate the motif environment
    conda activate motif
    
    echo "[APPEND_TO_DB] Starting database update process at \\\$(date)"
    echo "[APPEND_TO_DB] STREME results directory: ${streme_results}"
    echo "[APPEND_TO_DB] TRF results directory: ${trf_results}"
    
    # Check if a database file already exists in the project root directory
    DB_PATH="./motifs.db"
    if [ -f "\$DB_PATH" ]; then
        echo "[APPEND_TO_DB] Using existing database at \$DB_PATH"
        cp "\$DB_PATH" ./motifs.db
    else
        echo "[APPEND_TO_DB] No existing database found. Creating a new one."
        # Create a temporary script to create the database
        cat > create_db.sh << 'EOF'
#!/bin/bash

# Create a SQLite database to store the results of the motif discovery pipeline

DB_FILE="motifs.db"

# Check if the database file already exists
if [ -f "$DB_FILE" ]; then
    echo "Database file $DB_FILE already exists. Skipping creation."
    exit 0
fi

# Create the database and define the schema
echo "Creating new database file: $DB_FILE"
sqlite3 "$DB_FILE" <<EOFINNER
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
);

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
);
EOFINNER

echo "Database $DB_FILE created successfully with tables 'motifs' and 'metadata'."
EOF
        chmod +x create_db.sh
        ./create_db.sh
    fi
    
    # Check if database file exists now
    if [ -f "./motifs.db" ]; then
        echo "[APPEND_TO_DB] Database file found at ./motifs.db"
    else
        echo "[APPEND_TO_DB] ERROR: Database file not found at ./motifs.db"
        exit 1
    fi
    
    # Create a Python script to update the motif database
    echo "[APPEND_TO_DB] Creating Python script to update motif database"
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
    print("[APPEND_TO_DB] New database created at", db_path)

def update_motif_db(streme_dir, trf_dir, db_path):
    '''Update the motif database with STREME and TRF results'''
    # Create database if it doesn't exist
    if not os.path.exists(db_path):
        create_motif_db(db_path)
    
    conn = sqlite3.connect(db_path)
    print("[APPEND_TO_DB] Connected to database at", db_path)
    
    # Load STREME motifs
    streme_tsv = glob.glob(f"{streme_dir}/parsed_streme.tsv")
    if streme_tsv and os.path.exists(streme_tsv[0]) and os.path.getsize(streme_tsv[0]) > 0:
        try:
            streme_df = pd.read_csv(streme_tsv[0], sep='\t')
            streme_motif_count = len(streme_df)
            print(f"[APPEND_TO_DB] Loaded {streme_motif_count} motifs from STREME at {streme_tsv[0]}")
            
            # Get species name from directory structure
            species = os.path.basename(streme_dir).split('_')[0]
            streme_df['species'] = species
            print(f"[APPEND_TO_DB] Species for STREME motifs: {species}")
            
            # Insert STREME motifs into database
            streme_df.to_sql('motifs', conn, if_exists='append', index=False)
            print(f"[APPEND_TO_DB] Inserted {streme_motif_count} STREME motifs into database")
        except Exception as e:
            print(f"[APPEND_TO_DB] ERROR: Failed to load or insert STREME motifs: {e}")
            streme_motif_count = 0
    else:
        streme_motif_count = 0
        print("[APPEND_TO_DB] No STREME motifs to load")
        if streme_tsv:
            print(f"[APPEND_TO_DB] Checked file: {streme_tsv[0]}")
        else:
            print(f"[APPEND_TO_DB] No STREME TSV file found in {streme_dir}")
    
    # Load TRF motifs
    trf_tsv = glob.glob(f"{trf_dir}/parsed_trf.tsv")
    if trf_tsv and os.path.exists(trf_tsv[0]) and os.path.getsize(trf_tsv[0]) > 0:
        try:
            trf_df = pd.read_csv(trf_tsv[0], sep='\t')
            trf_motif_count = len(trf_df)
            print(f"[APPEND_TO_DB] Loaded {trf_motif_count} motifs from TRF at {trf_tsv[0]}")
            
            # Get species name from directory structure
            species = os.path.basename(trf_dir).split('_')[0]
            trf_df['species'] = species
            print(f"[APPEND_TO_DB] Species for TRF motifs: {species}")
            
            # Insert TRF motifs into database
            trf_df.to_sql('motifs', conn, if_exists='append', index=False)
            print(f"[APPEND_TO_DB] Inserted {trf_motif_count} TRF motifs into database")
        except Exception as e:
            print(f"[APPEND_TO_DB] ERROR: Failed to load or insert TRF motifs: {e}")
            trf_motif_count = 0
    else:
        trf_motif_count = 0
        print("[APPEND_TO_DB] No TRF motifs to load")
        if trf_tsv:
            print(f"[APPEND_TO_DB] Checked file: {trf_tsv[0]}")
        else:
            print(f"[APPEND_TO_DB] No TRF TSV file found in {trf_dir}")
    
    # Update metadata
    import datetime
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[APPEND_TO_DB] Updating metadata with timestamp: {now}")
    
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
        print(f"[APPEND_TO_DB] Window stats - Total: {window_count}, Filtered: {filtered_window_count}")
    else:
        print(f"[APPEND_TO_DB] No window stats file found")
    
    # Get species name from one of the result directories
    species = os.path.basename(streme_dir).split('_')[0] or os.path.basename(trf_dir).split('_')[0]
    print(f"[APPEND_TO_DB] Species for metadata: {species}")
    
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
        print(f"[APPEND_TO_DB] Updated metadata for species {species}")
    else:
        conn.execute('''
        INSERT INTO metadata
        (species, window_count, filtered_window_count, streme_motif_count, trf_motif_count, processing_date)
        VALUES (?, ?, ?, ?, ?, ?)
        ''', (species, window_count, filtered_window_count, streme_motif_count, trf_motif_count, now))
        print(f"[APPEND_TO_DB] Inserted new metadata for species {species}")
    
    conn.commit()
    conn.close()
    print(f"[APPEND_TO_DB] Updated motif database with {streme_motif_count} STREME motifs and {trf_motif_count} TRF motifs")
    print(f"[APPEND_TO_DB] Database connection closed")

if __name__ == "__main__":
    streme_dir = sys.argv[1]
    trf_dir = sys.argv[2]
    db_path = sys.argv[3]
    
    update_motif_db(streme_dir, trf_dir, db_path)
EOF
    
    # Make the script executable
    echo "[APPEND_TO_DB] Making update_motif_db.py script executable"
    chmod +x update_motif_db.py
    
    # Run the database update script
    echo "[APPEND_TO_DB] Running database update script"
    python update_motif_db.py ${streme_results} ${trf_results} motifs.db
    if [ \$? -eq 0 ]; then
        echo "[APPEND_TO_DB] Database update script completed successfully"
    else
        echo "[APPEND_TO_DB] ERROR: Database update script failed"
    fi
    echo "[APPEND_TO_DB] Database update process completed on \\\$(date)"
    """
}