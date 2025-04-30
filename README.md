# motif_graph_explore



steps:



PrePipeline:
1. Scrape the ncbi data
    script: ncbi_scrape.py

Pipeline:
1. Download genome
    script: download_genome.nf
2. Index genome 
    script: index_genome.nf
3. softmask_genome
    script: softmask_genome.nf
4. make and extract_windows
    script: make_extract_windows.nf

4. run STREME analysis
    script: STREME.nf
4. run TRF 
    script: TRF.nf
5. append results to db
    script: append_to_db.nf
6. clean_temp
    script: clean_temp.nf

PostPipeline:
