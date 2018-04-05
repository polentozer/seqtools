#!/bin/bash
echo -e "\nTesting...\n"
mkdir ./temp_test
touch ./log
test=0
python3 seqtools -t test_data/table.csv -i test_data/prot.fasta -o temp_test/prot.fasta -pf
python3 seqtools -t test_data/table.csv -i test_data/multi_prot.fasta -o temp_test/multi_prot.fasta -pf
python3 seqtools -t test_data/table.csv -i test_data/dna.fasta -o temp_test/dna.fasta -f
python3 seqtools -t test_data/table.csv -i test_data/multi_dna.fasta -o temp_test/multi_dna.fasta -f
for file in ./temp_test/*
do
    if cmp -s "${file}" "./test_data/solution_${file:12}"; then 
        echo -e "${file:12}: OK" >> ./log
    else
        test=1
        echo -e "${file:12}: BAD" >> ./log
    fi
done
if [ $test -eq 0 ]; then
    echo -e "\nAll's well."
    rm -r temp_test log
else
    echo -e "\nOops, there were some problems."
fi
