#!/usr/bin/env python3
# Mandatory imports
import argparse
import os

try:
    import pandas
except ImportError:
    raise ImportError("Missing package: `PANDAS`!")


def reverse_translate(amino, codon_table):
    """
    Reverse-translates aminoacid back to the most used codon from codon table.
    """
    temp_values = []
    temp_triplets = []

    for i, row in codon_table.loc[codon_table[0] == amino].iterrows():
        temp_triplets.append((row[1]).strip().upper())
        temp_values.append(row[2])
    
    return temp_triplets[temp_values.index(max(temp_values))]


def make_dna_triplets(string):
    """
    Makes chunks of 3 characters from a long string.
    """
    return [string[start:start+3] for start in range(0, len(string), 3)]


def codon_optimize(name, orf, codon_table):
    """
    Checks if a given sequence starts with the start codon ('ATG') and promts you if it doesn't.
    Optimizes codons to most frequently used from codon table.
    """
    forced = False

    if orf[:3] != "ATG":
        cont_prompt = input("Sequence with id `{0}` is not ORF, do you want to continue anyway? (Y/n): ".format(name))
    
        if cont_prompt.upper() == "Y" or cont_prompt.upper() == "YES" or cont_prompt == "":
            forced = True
        else:
            return False, forced

    optimized_orf = ""

    for triplet in make_dna_triplets(orf):
        
        if len(triplet) == 3:
            optimized_orf += reverse_translate(codon_table.loc[
                codon_table[1] == triplet][0].iloc[0][0], codon_table)
        else:
            optimized_orf += triplet

    return optimized_orf, forced


def writer(string, path):
    """
    Simple function for writing `.fasta` files.
    """
    with open(path, 'w') as file_out:
        file_out.write(string)


def open_fasta(file_input_path):
    """
    Function for opening `.fasta` files that can contain more than just one sequence.
    Returns dictionary of sequences ==> {id: sequence}
    """
    sequence_dictionary = {}

    with open(file_input_path, 'r') as input_file:
        data = input_file.readlines()
    
    for line in data:
        temp = line.strip()
        if line[0] == ">":
            sequence_id = temp[:20]
            sequence_dictionary[sequence_id] = ""
        else:
            sequence_dictionary[sequence_id] += temp.upper()

    return sequence_dictionary


def main():
    """
    Main function with argument parsing.
    """
    description = ("Welcome to python script for converting protein sequence to DNA sequence with\
                    custom table for codon usage.")

    # Arguments
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-p","--protein", action="store_true", help="Use this flag when working with protein sequences", required=False)
    parser.add_argument("-f","--file", help="Path to target fasta file", required=True)
    parser.add_argument("-t","--table",help="Path to codon usage table in csv format: 'aminoacid,triplet,value')", required=True)
    parser.add_argument("-s", "--save", help="Path for output fasta file", required=False)
    args = parser.parse_args()

    # Variables
    current_dir = os.path.dirname(os.path.realpath(__file__))
    path_fasta = "{0}/{1}".format(current_dir, args.file)
    path_table = "{0}/{1}".format(current_dir, args.table)
    codon_table = pandas.read_csv(path_table, header=None)
    sequences = open_fasta(path_fasta)

    if args.save:
        path_save = "{0}/{1}".format(current_dir, args.save)
        msg_saved = "Output saved to `{0}`".format(args.save)

    if args.protein:
        reverse_translations = {}

        for id, sequence in sequences.items():
            reverse_translations[id] = ""
            for amino in sequence:
                reverse_translations[id] += reverse_translate(amino, codon_table)

        string = "\n".join(["{0}|OPTIMIZED\n{1}\n".format(k, v) for k, v in reverse_translations.items()])

        if args.save:
            writer(string, path_save)
            print(msg_saved)
        else:
            print("\nCodon optimized reverse translation:\n\n{0}".format(string))
    
    else:
        optimized_orfs = {}

        for id, sequence in sequences.items():
            optimized_orfs[id] = codon_optimize(id, sequence, codon_table)

        string = ""
        for k, v in optimized_orfs.items():
            if v[0]:
                k += "|OPTIMIZED"
                if v[1]:
                    k += "|forced"
                string += "{0}\n{1}\n\n".format(k, v[0])

        if args.save:
            writer(string, path_save)
            print(msg_saved)
        else:
            print("\nCodon optimized ORFs:\n\n{0}".format(string))


if __name__ == "__main__":
    main()