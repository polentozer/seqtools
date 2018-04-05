#!/bin/sh
echo "testing"
mkdir ./temp_test
seqtools -p -t test_data/table.csv -f test_data/dna.fasta -s temp_test/temp_dna.fasta
seqtools -p -t test_data/table.csv -f test_data/dna.fasta -s temp_test/temp_dna.fasta
seqtools -t test_data/table.csv -f test_data/dna.fasta -s temp_test/temp_dna.fasta
seqtools -t test_data/table.csv -f test_data/multi_dna.fasta -s temp_test/temp_multi_dna.fasta