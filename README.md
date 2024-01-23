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