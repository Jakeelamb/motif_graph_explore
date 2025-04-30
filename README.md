# motif_graph_explore



clone the repo
```
git clone 


create the motif conda environment
```
conda env create -f environment.yml
```










steps:

PrePipeline:
1. Scrape the ncbi data
    script: ncbi_scrape.py

Pipeline:
1. Download genome
    script: download_genome.nf
2. Index genome (depends on download_genome)
    script: index_genome.nf
3. softmask_genome (depends on index_genome)
    script: softmask_genome.nf
4. make and extract_windows (depends on softmask_genome)
    script: make_extract_windows.nf
5. run STREME analysis (depends on make_extract_windows)
    script: STREME.nf
3. run TRF (depends on index_genome)
    script: TRF.nf
5. append results to db (depends on TRF and STEME)
    script: append_to_db.nf
6. clean_temp (depends on append_to_db)
    script: clean_temp.nf

PostPipeline:
