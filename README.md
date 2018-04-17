# SEQTOOLS
## Description
Command line tool for manipulating `.fasta` files.

## Installation
Script can be run on it's own or it could be installed with python's `setup.py` script.

### Requirements
- pandas

Codon usage table should be provided in a `.csv` file in the following format:
```csv
aminoacid,dna_triplet,value
aminoacid,dna_triplet,value
aminoacid,dna_triplet,value
...
```

To install type in terminal:
```bash
python setup.py install
```

## Usage
```bash
# translate DNA sequence to protein sequence using genetic code table provided in the codon usage table and print the output
seqtools -t test_data/table.csv -i test_data/dna.fasta

# reverse-translate and codon-optimize a protein sequence and print the output
seqtools -t test_data/table.csv -i test_data/prot.fasta -p

# codon-optimize a DNA sequence and save the output to a file
seqtools -O -t test_data/table.csv -i test_data/dna.fasta  -o test.fasta

# force codon-optimization on a multi-DNA-fasta file with sequences that don't
# start with the 'ATG' triplet (>name|OPTIMIZED|forced)
seqtools -Oft test_data/table.csv -i test_data/multi_dna.fasta
```