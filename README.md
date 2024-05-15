# TE_Idx

## Commands
### bgzf-filter 
- --assembly <ASSEMBLY> : Name of assembly/assembly folder
- --data-type <DATA-TYPE> : Type of data to be searched \(Assembly, Benchmarks, Masks)
- --fam <FAM> : Family name, corresponds to compressed TSV file prefix
- --position <POSITION> : number corresponding to the search field (column), 1-indexed
- --term <TERM> : Term to be searched for. If absent, all rows will be returned
- --outfile <OUTFILE> : filename for output bgzf compressed TSV output file
- --web-fmt <WEB-FMT> : Flag to reformat the feilds to match Dfam.org download file format

### build-idx 
- --assembly <ASSEMBLY> : path to assembly folder to index
- --data-type <DATA-TYPE> : type of data to index \(Assembly, Benchmarks, Masks)

### prep-beds
- --assembly <ASSEMBLY> : path to assembly folder to index
- --in-tsv <IN-TSV> : input file
- --data-type <DATA-TYPE> : type of data to index \(Assembly, Benchmarks, Masks)

### process-json
- --in-file <> : Input JSON file, exported from PHPMyAdmin
- --key <> : Key attribute, such as accession, to be used in the new map
- --outfile <> : Optional: Output file. If not provided, will print to stdout

### prepare-assembly
- --assembly <ASSEMBLY> : path to assembly folder to index

### idx-query
- --assembly <ASSEMBLY> : path to assembly folder to index
- --data-type <DATA-TYPE> : type of data to index \(Assembly, Benchmarks, Masks)
- --chrom <CHROM> : chromosome number/accession
- --start <START> : start position
- --end <END> : end position
- --family <FAMILY> : Optional: Only return hits matching accession
- --nrph <NRPH> : Only return NRPH hits

### json-query
- --assembly <ASSEMBLY> : path to assembly folder to index
- --data-type <DATA-TYPE> : type of data to index \(Assembly, Benchmarks, Masks)
- --key <KEY> : Key value to search by
- --target <TARGET> : Target Attribute to return

### read-family-assembly-annotations
- --id <ID> : Family Accession
- --assembly-id <ASSEMBLY-ID> : Name of assembly/assembly folder
- --nrph <NRPH> : Only Return NRPH hits
- --outfile <OUTFILE> : Optional: Output file

# Sources
* hg38-byacc-bench_region.tsv -> buildFullRegion.py.test
* hg38-byacc-full_region.tsv -> buildFullRegion.py.test
* hg38-byseq-bench_region.tsv -> buildFullRegion.py.test
* hg38-byseq-full_region.tsv -> buildFullRegion.py.test
* hg38-mask.tsv -> SELECT * FROM `mask`
* hg38-model_lengths.json -> SELECT t2.family_accession, t1.length FROM dfam_dev.family as t1 JOIN hg38_df38.model_file as t2 ON t1.accession=t2.family_accession
* hg38-sequence.json -> SELECT * FROM `sequence`

## Testing
`cargo test`