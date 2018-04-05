## SEQTOOLS
### Description
Command line tool for manipulating `.fasta` files.

#### Requirements
pandas

table.csv has to be in the following style for now:
```csv
aminoacid,dna_triplet,value
aminoacid,dna_triplet,value
aminoacid,dna_triplet,value
...
```

### Usage
```bash
# reverse-translate and codon-optimize a protein sequence and print the output
./seqtools -pf test_data/prot.fasta -t test_data/table.csv

# codon-optimize a DNA sequence and save the output to a file
./seqtools -f test_data/dna.fasta -t test_data/table.csv -s test.fasta
```