// run STREME analysis on the 10kb fasta windows 

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
    source /nfs/home/jlamb/bin/miniconda3/etc/profile.d/conda.sh
    
    # Activate the motif environment
    conda activate motif
    
    echo "[STREME] Starting STREME analysis for: ${windows} at \\\$(date)"
    echo "[STREME] DEBUGGING MODE: Using RNA alphabet for faster computation"
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
    
    # Run STREME with optimized parameters for better performance
    # Using RNA mode and fixed seed for reproducibility during debugging
    echo "[STREME] Running STREME analysis with optimized parameters"
    streme --p ${windows}/filtered_windows.fa \
        --rna \
        --oc ${windows}_streme/streme_out \
        --minw 6 \
        --maxw 12 \
        --order 2 \
        --thresh 0.01 \
        --neval 15 \
        --nref 2 \
        --niter 10 \
        --nmotifs 75 \
        --hofract 0.05 \
        --totallength 10000000 \
        --time 86400 \
        --seed 42 \
        --no-pgc
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
