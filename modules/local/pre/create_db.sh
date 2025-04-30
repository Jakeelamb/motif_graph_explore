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
sqlite3 "$DB_FILE" <<EOF
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
EOF

echo "Database $DB_FILE created successfully with tables 'motifs' and 'metadata'."