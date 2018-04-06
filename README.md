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
./seqtools -t test_data/table.csv -i test_data/prot.fasta -p

# codon-optimize a DNA sequence and save the output to a file
./seqtools -t test_data/table.csv -i test_data/dna.fasta  -o test.fasta

# force codon-optimization on a multi-DNA-fasta file with sequences that don't
# start with the 'ATG' triplet (>name|OPTIMIZED|forced)
./seqtools -t test_data/table.csv -i test_data/multi_dna.fasta -f
```