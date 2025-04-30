#!/usr/bin/env nextflow

/*
=============================
MOTIF DISCOVERY PIPELINE
=============================
*/

// Input parameters
params.input = 'accensions.txt'
params.outdir = './results'

// Print the pipeline header
log.info """\
         MOTIF DISCOVERY PIPELINE
         ===================================
         Input file         : ${params.input}
         Output directory   : ${params.outdir}
         """
         .stripIndent()

// Create the channel for accession numbers
Channel
    .fromPath(params.input)
    .splitCsv(header:true, sep:'\t')
    .map { row -> row.accension }
    .set { accessions_ch }

// Define the main workflow
workflow {
    // Pre-processing
    // Empty for now - could include ncbi_scrape.py

    // Pipeline execution
    // Step 1: Download genome
    fetch_genome(accessions_ch)
    
    // Step 2: Index genome
    index_genome(fetch_genome.out)
    
    // Step 3: Softmask genome
    softmask_genome(index_genome.out)
    
    // Step 4: Make and extract windows
    make_extract_windows(softmask_genome.out)
    
    // Step 5: Run STREME analysis
    STREME(make_extract_windows.out)
    
    // Step 6: Run TRF
    TRF(index_genome.out)
    
    // Step 7: Append results to database
    append_to_db(STREME.out, TRF.out)
    
    // Step 8: Clean temporary files
    clean_temp(append_to_db.out)
    
    // Post-processing
    // Empty for now
}

// Include the processes
include { fetch_genome } from './modules/local/pipeline/fetch_genome'
include { index_genome } from './modules/local/pipeline/index_genome'
include { softmask_genome } from './modules/local/pipeline/softmask_genome'
include { make_extract_windows } from './modules/local/pipeline/make_extract_windows'
include { STREME } from './modules/local/pipeline/STREME'
include { TRF } from './modules/local/pipeline/TRF'
include { append_to_db } from './modules/local/pipeline/append_to_db'
include { clean_temp } from './modules/local/pipeline/clean_temp'
