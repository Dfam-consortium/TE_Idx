# TE_Idx

## Commands
### bgzf-filter 
--infile <INFILE> : bgzf compressed TSV (whitespace separated) input file
--position <POSITION> : number corresponding to the search field (column)
--term <TERM> : search term to be matched
--outfile <OUTFILE> : filename for output bgzf compressed TSV output file

### find-family 
--id <ID> 
--assembly <ASSEMBLY>

### idx 
--assembly <ASSEMBLY> : path to assembly folder to index
--build <BUILD> : boolean, triggers the construction of an index
--search <SEARCH> : String in the form of Chromosome-Start-End coordinates to search

# Sources
* hg38-byacc-bench_region.tsv -> buildFullRegion.py.test
* hg38-byacc-full_region.tsv -> buildFullRegion.py.test
* hg38-byseq-bench_region.tsv -> buildFullRegion.py.test
* hg38-byseq-full_region.tsv -> buildFullRegion.py.test
* hg38-mask.tsv -> SELECT * FROM `mask`
* hg38-model_lengths.json -> SELECT t2.family_accession, t1.length FROM dfam_dev.family as t1 JOIN hg38_df38.model_file as t2 ON t1.accession=t2.family_accession
* hg38-sequence.json -> SELECT * FROM `sequence`