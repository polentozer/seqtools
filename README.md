# SEQTOOLS
## Description
Command line tool for manipulating `.fasta` files.

## Installation
Tool can be installed by building it first with `python setup.py sdist` and then installing it with `pip install .`.

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
# translate DNA sequence to protein sequence using genetic code table provided in the default table and print the output
seqtools -i test_data/sample_dna.fasta

# reverse-translate and codon-optimize a protein sequence and print the output
seqtools -t test_data/sample_table.csv -i test_data/sample_prot.fasta -p

# codon-optimize a DNA sequence and save the output to a file
seqtools -O -t test_data/sample_table.csv -i test_data/sample_dna.fasta  -o test.fasta

# force codon-optimization on a multi-DNA-fasta file with sequences that don't
# start with the 'ATG' triplet (>name|OPTIMIZED|forced)
seqtools -Oft test_data/sample_table.csv -i test_data/sample_multi_dna.fasta
```