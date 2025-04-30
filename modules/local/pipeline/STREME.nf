// run STREME analysis on the 10kb fasta windows 

// main.nf - Motif Discovery Pipeline

workflow {
  fetch_genome()
  mask_genome(fetch_genome.out)
  index_genome(mask_genome.out)
  make_windows(index_genome.out)
  extract_windows(make_windows.out)
  run_trf(mask_genome.out)
  run_streme(extract_windows.out)
  parse_streme(run_streme.out)
  parse_trf(run_trf.out)
  update_motif_db(parse_streme.out, parse_trf.out)
}

process fetch_genome {
  input:
  val accession

  output:
  path "genome/${accession}.fa.gz" into fetch_genome_out

  script:
  """
  conda activate motif
  datasets download genome accession $accession --filename ${accession}.zip
  unzip ${accession}.zip -d genome/
  """
}

process mask_genome {
  input:
  path genome_file from fetch_genome_out

  output:
  path "${genome_file.baseName}.masked.fa" into mask_genome_out

  script:
  """
  conda activate motif
  RepeatMasker -xsmall -pa 4 -species vertebrate -dir . ${genome_file}
  mv ${genome_file.baseName}.masked ${genome_file.baseName}.masked.fa
  """
}

process index_genome {
  input:
  path genome from mask_genome_out

  output:
  path "${genome}.fai" into index_genome_out

  script:
  """
  conda activate motif
  samtools faidx ${genome}
  """
}

process make_windows {
  input:
  path fai from index_genome_out

  output:
  path "windows.bed" into make_windows_out

  script:
  """
  conda activate motif
  bedtools makewindows -g ${fai} -w 10000 > windows.bed
  """
}

process extract_windows {
  input:
  path bed from make_windows_out
  path fasta from mask_genome_out

  output:
  path "genome.windows.fa" into extract_windows_out

  script:
  """
  conda activate motif
  bedtools getfasta -fi ${fasta} -bed ${bed} -fo genome.windows.fa
  """
}

process run_streme {
  input:
  path windows_fa from extract_windows_out

  output:
  path "streme_out" into run_streme_out

  script:
  """
  # Initialize Conda
  source /home/jake/miniconda3/etc/profile.d/conda.sh
  
  # Activate the motif environment
  conda activate motif
  
  echo "[STREME] Starting STREME analysis for: ${windows} at \\\$(date)"
  mkdir -p ${windows}_streme
  
  # Check if input file exists and has data
  if [ ! -f "${windows}/filtered_windows.fa" ]; then
      echo "[STREME] ERROR: Input FASTA file not found at ${windows}/filtered_windows.fa"
      exit 1
  else
      echo "[STREME] Input FASTA file found at ${windows}/filtered_windows.fa"
      echo "[STREME] Number of sequences in input: \$(grep -c '>' ${windows}/filtered_windows.fa)"
      echo "[STREME] File size: \$(ls -lh ${windows}/filtered_windows.fa | awk '{print \$5}')"
  fi
  
  # Run STREME on the full set of windows without subsampling
  echo "[STREME] Running STREME analysis on full dataset"
  streme --p ${windows}/filtered_windows.fa \
      --dna \
      --oc ${windows}_streme/streme_out \
      --minw 5 \
      --maxw 15 \
      --order 2 \
      --thresh 0.01 \
      --nmotifs 100 \
      --no-pgc \
      --totallength 1000000000
  if [ \$? -eq 0 ]; then
      echo "[STREME] STREME analysis completed successfully"
  else
      echo "[STREME] ERROR: STREME analysis failed"
  fi
  
  # Check if STREME output file exists
  if [ -f "${windows}_streme/streme_out/streme.txt" ]; then
      echo "[STREME] STREME output file found at ${windows}_streme/streme_out/streme.txt"
      echo "[STREME] Number of motifs in output: \$(grep -c '^MOTIF' ${windows}_streme/streme_out/streme.txt)"
  else
      echo "[STREME] ERROR: STREME output file not found at ${windows}_streme/streme_out/streme.txt"
  fi
      
  # Parse STREME results to TSV format
  echo "[STREME] Parsing STREME results to TSV format"
  python3 -c "import os; import sys; import re; import pandas as pd; streme_txt_path = '${windows}_streme/streme_out/streme.txt'; output_tsv = '${windows}_streme/parsed_streme.tsv'; motifs = []; with open(streme_txt_path, 'r') if os.path.exists(streme_txt_path) else open('/dev/null', 'r') as f: lines = f.readlines(); in_motif_section = False; current_motif = {}; for line in lines: line = line.strip(); if line.startswith('MOTIF'): in_motif_section = True; motif_match = re.match(r'MOTIF\\s+(\\d+)', line); if motif_match: current_motif = {'motif_id': f'STREME-{motif_match.group(1)}'}; elif in_motif_section and line.startswith('p-value'): pvalue_match = re.match(r'p-value\\s+=\\s+(\\S+)', line); if pvalue_match: current_motif['p_value'] = pvalue_match.group(1); elif in_motif_section and line.startswith('E-value'): evalue_match = re.match(r'E-value\\s+=\\s+(\\S+)', line); if evalue_match: current_motif['e_value'] = evalue_match.group(1); elif in_motif_section and line.startswith('Consensus'): consensus_match = re.match(r'Consensus\\s+=\\s+(\\S+)', line); if consensus_match: current_motif['consensus'] = consensus_match.group(1); elif in_motif_section and line.startswith('Width'): width_match = re.match(r'Width\\s+=\\s+(\\d+)', line); if width_match: current_motif['width'] = width_match.group(1); elif in_motif_section and line.startswith('Sites'): sites_match = re.match(r'Sites\\s+=\\s+(\\d+)', line); if sites_match: current_motif['sites'] = sites_match.group(1); motifs.append(current_motif); in_motif_section = False; current_motif = {}; if motifs: df = pd.DataFrame(motifs); df['source'] = 'STREME'; df.to_csv(output_tsv, sep='\\t', index=False); print(f'[STREME] Parsed {len(motifs)} motifs from STREME output'); else: print('[STREME] No motifs found in STREME output'); columns = ['motif_id', 'consensus', 'width', 'sites', 'p_value', 'e_value', 'source']; pd.DataFrame(columns=columns).to_csv(output_tsv, sep='\\t', index=False);"
  echo "[STREME] STREME process completed on \\\$(date)"
  """
}

process run_trf {
  input:
  path genome_fa from mask_genome_out

  output:
  path "trf.dat" into run_trf_out

  script:
  """
  conda activate motif
  trf-big ${genome_fa} 2 7 7 80 10 50 500 -d -h
  mv ${genome_fa}.2.7.7.80.10.50.500.dat trf.dat
  """
}

process parse_streme {
  input:
  path streme_dir from run_streme_out

  output:
  path "parsed_streme.tsv" into parse_streme_out

  script:
  """
  conda activate motif
  python3 scripts/parse_streme.py --input ${streme_dir}/streme.txt --output parsed_streme.tsv
  """
}

process parse_trf {
  input:
  path trf_dat from run_trf_out

  output:
  path "parsed_trf.tsv" into parse_trf_out

  script:
  """
  conda activate motif
  python3 scripts/parse_trf.py --input ${trf_dat} --output parsed_trf.tsv
  """
}

process update_motif_db {
  input:
  path streme_tsv from parse_streme_out
  path trf_tsv from parse_trf_out

  output:
  path "motif.db" emit: db_out

  script:
  """
  conda activate motif
  python3 scripts/update_motif_db.py --streme ${streme_tsv} --trf ${trf_tsv} --db motif.db
  """
}

process STREME {
    tag "$windows"
    publishDir "${params.outdir}/streme", mode: 'copy'
    
    input:
    path windows
    
    output:
    path "${windows}_streme", emit: streme_results
    
    script:
    """
    # Initialize Conda
    source /home/jake/miniconda3/etc/profile.d/conda.sh
    
    # Activate the motif environment
    conda activate motif
    
    echo "[STREME] Starting STREME analysis for: ${windows} at \\\$(date)"
    mkdir -p ${windows}_streme
    
    # Check if input file exists and has data
    if [ ! -f "${windows}/filtered_windows.fa" ]; then
        echo "[STREME] ERROR: Input FASTA file not found at ${windows}/filtered_windows.fa"
        exit 1
    else
        echo "[STREME] Input FASTA file found at ${windows}/filtered_windows.fa"
        echo "[STREME] Number of sequences in input: \$(grep -c '>' ${windows}/filtered_windows.fa)"
        echo "[STREME] File size: \$(ls -lh ${windows}/filtered_windows.fa | awk '{print \$5}')"
    fi
    
    # Run STREME on the full set of windows without subsampling
    echo "[STREME] Running STREME analysis on full dataset"
    streme --p ${windows}/filtered_windows.fa \
        --dna \
        --oc ${windows}_streme/streme_out \
        --minw 5 \
        --maxw 15 \
        --order 2 \
        --thresh 0.01 \
        --nmotifs 100 \
        --no-pgc \
        --totallength 1000000000
    if [ \$? -eq 0 ]; then
        echo "[STREME] STREME analysis completed successfully"
    else
        echo "[STREME] ERROR: STREME analysis failed"
    fi
    
    # Check if STREME output file exists
    if [ -f "${windows}_streme/streme_out/streme.txt" ]; then
        echo "[STREME] STREME output file found at ${windows}_streme/streme_out/streme.txt"
        echo "[STREME] Number of motifs in output: \$(grep -c '^MOTIF' ${windows}_streme/streme_out/streme.txt)"
    else
        echo "[STREME] ERROR: STREME output file not found at ${windows}_streme/streme_out/streme.txt"
    fi
        
    # Parse STREME results to TSV format
    echo "[STREME] Parsing STREME results to TSV format"
    python3 -c "import os; import sys; import re; import pandas as pd; streme_txt_path = '${windows}_streme/streme_out/streme.txt'; output_tsv = '${windows}_streme/parsed_streme.tsv'; motifs = []; with open(streme_txt_path, 'r') if os.path.exists(streme_txt_path) else open('/dev/null', 'r') as f: lines = f.readlines(); in_motif_section = False; current_motif = {}; for line in lines: line = line.strip(); if line.startswith('MOTIF'): in_motif_section = True; motif_match = re.match(r'MOTIF\\s+(\\d+)', line); if motif_match: current_motif = {'motif_id': f'STREME-{motif_match.group(1)}'}; elif in_motif_section and line.startswith('p-value'): pvalue_match = re.match(r'p-value\\s+=\\s+(\\S+)', line); if pvalue_match: current_motif['p_value'] = pvalue_match.group(1); elif in_motif_section and line.startswith('E-value'): evalue_match = re.match(r'E-value\\s+=\\s+(\\S+)', line); if evalue_match: current_motif['e_value'] = evalue_match.group(1); elif in_motif_section and line.startswith('Consensus'): consensus_match = re.match(r'Consensus\\s+=\\s+(\\S+)', line); if consensus_match: current_motif['consensus'] = consensus_match.group(1); elif in_motif_section and line.startswith('Width'): width_match = re.match(r'Width\\s+=\\s+(\\d+)', line); if width_match: current_motif['width'] = width_match.group(1); elif in_motif_section and line.startswith('Sites'): sites_match = re.match(r'Sites\\s+=\\s+(\\d+)', line); if sites_match: current_motif['sites'] = sites_match.group(1); motifs.append(current_motif); in_motif_section = False; current_motif = {}; if motifs: df = pd.DataFrame(motifs); df['source'] = 'STREME'; df.to_csv(output_tsv, sep='\\t', index=False); print(f'[STREME] Parsed {len(motifs)} motifs from STREME output'); else: print('[STREME] No motifs found in STREME output'); columns = ['motif_id', 'consensus', 'width', 'sites', 'p_value', 'e_value', 'source']; pd.DataFrame(columns=columns).to_csv(output_tsv, sep='\\t', index=False);"
    echo "[STREME] STREME process completed on \\\$(date)"
    """
}
