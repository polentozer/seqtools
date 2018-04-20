import argparse
import os
import pandas
import seqtools

def main():
    """
    Main function with argument parsing.
    """
    description = ("Welcome to python script for converting protein sequence to DNA sequence with\
                    custom table for codon usage.")

    # Arguments
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-v", "--version", action="version", version=seqtools.__version__)
    parser.add_argument("-i", "--input", help="Path to input 'fasta' file", required=True)
    parser.add_argument("-t", "--table",help="Path to codon usage table in csv format: 'aminoacid,triplet,value')", required=True)
    parser.add_argument("-p", "--protein", action="store_true", help="Use this flag when working with protein sequences", required=False)
    parser.add_argument("-o",  "--output", help="Name for the output fasta file", required=False)
    parser.add_argument("-f", "--force", action="store_true", help="Use this flag to omit any prompts", required=False)
    parser.add_argument("-O", "--optimize", action="store_true", help="Use to optimize DNA sequence instead translating it.", required=False)
    args = parser.parse_args()

    # Variables
    path_fasta = os.path.join(os.path.join(os.getcwd(), args.input))
    path_table = os.path.join(os.path.join(os.getcwd(), args.table))
    codon_table = pandas.read_csv(path_table, header=None)
    sequences = seqtools.io.fasta.open_fasta(path_fasta)

    # Protein translation
    if args.protein:
        solution = seqtools.seqtools.protein_to_dna(sequences, codon_table)

    elif not args.protein:
        solution = seqtools.seqtools.dna_operation(sequences, codon_table, args.force, args.optimize)

    # Saving/Printing solution(s)
    if args.output:
        path_save = os.path.join(os.path.join(os.getcwd(), args.output))
        msg_saved = "Output saved to `{0}`".format(args.output)
        seqtools.io.output.writer(solution, path_save)
        print(msg_saved)
    else:
        print()
        print(seqtools.io.output.writer(solution))

