#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

cd "$parent_path"

echo -e "\nTesting...\n"

mkdir ./temp_test
touch ./log
test=0

seqtools -t ../test_data/sample_table.csv -i ../test_data/sample_prot.fasta -o temp_test/prot.fasta -pfO
seqtools -t ../test_data/sample_table.csv -i ../test_data/sample_multi_prot.fasta -o temp_test/multi_prot.fasta -pfO
seqtools -t ../test_data/sample_table.csv -i ../test_data/sample_dna.fasta -o temp_test/dna.fasta -fO
seqtools -t ../test_data/sample_table.csv -i ../test_data/sample_multi_dna.fasta -o temp_test/multi_dna.fasta -fO

for file in ./temp_test/*
do
    if cmp -s "${file}" "../test_data/solution_${file:12}"; then 
        echo -e "${file:12}: OK" >> ./log
    else
        test=1
        echo -e "${file:12}: BAD" >> ./log
    fi
done

echo -e "OVER\n" >> ./log

if [ $test -eq 0 ]; then
    echo -e "\nAll's well."
    rm -r temp_test log
else
    echo -e "\nOops, there were some problems."
fi
