# get the file size of the genome assembly and the taxonomy information from the metazoan_chromosome_assemblies_with_lineage.csv

import pandas as pd
import os
import sqlite3

csv_file = 'metazoan_chromosome_assemblies_with_lineage.csv'

if not os.path.exists(csv_file):
    print(f"Error: The file {csv_file} was not found. Please ensure it is in the current directory.")
    exit(1)

# Read the CSV file
df = pd.read_csv(csv_file)

# Connect to the database
conn = sqlite3.connect('motifs.db')
cur = conn.cursor()

# Create a table for taxonomy and size if it doesn't exist
cur.execute('''
    CREATE TABLE IF NOT EXISTS species (
        species_id INTEGER PRIMARY KEY AUTOINCREMENT,
        organism TEXT UNIQUE,
        taxonomy TEXT,
        genome_size_mb REAL
    )
''')

# Add new columns if they don't exist
for column in ['kingdom', 'phylum', 'class', 'taxonomic_order', 'family', 'genus', 'species_name', 'species_taxid', 'accession']:
    try:
        cur.execute(f'ALTER TABLE species ADD COLUMN {column} TEXT')
    except sqlite3.OperationalError:
        pass  # Column already exists

# Iterate through the CSV rows and update the database
for index, row in df.iterrows():
    organism = row['organism']
    kingdom = row.get('kingdom', '')
    phylum = row.get('phylum', '')
    class_ = row.get('class', '')
    taxonomic_order = row.get('order', '')
    family = row.get('family', '')
    genus = row.get('genus', '')
    species_name = row.get('species', '')
    species_taxid = row.get('speciestaxid', '')
    taxonomy = row.get('lineage', '')
    file_path = row.get('file_path', '')
    accession = row.get('accession', '')
    
    if os.path.exists(file_path):
        size_mb = os.path.getsize(file_path) / (1024 * 1024)  # Convert bytes to MB
    else:
        size_mb = 0.0
    
    cur.execute('''
        INSERT OR REPLACE INTO species (
            organism, kingdom, phylum, class, taxonomic_order, family, genus, 
            species_name, species_taxid, taxonomy, genome_size_mb, accession
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', (organism, kingdom, phylum, class_, taxonomic_order, family, genus, species_name, species_taxid, taxonomy, size_mb, accession))

conn.commit()
conn.close()

print('Taxonomy and genome size information updated in motifs.db.')

