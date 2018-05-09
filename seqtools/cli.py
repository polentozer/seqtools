import argparse
from os import path, getcwd
import pandas
import seqtools

def main():
    """
    Main function with argument parsing.
    """

    # Arguments
    parser = argparse.ArgumentParser(description=seqtools.__summary__)
    parser.add_argument("-v", "--version", action="version", version=seqtools.__version__)
    parser.add_argument("-i", "--input", help="Path to input 'fasta' file", required=True)
    parser.add_argument("-t", "--table",help="Path to codon usage table in csv format: 'aminoacid,triplet,value')", required=False)
    parser.add_argument("-p", "--protein", action="store_true", help="Use this flag when working with protein sequences", required=False)
    parser.add_argument("-o",  "--output", help="Name for the output fasta file", required=False)
    parser.add_argument("-f", "--force", action="store_true", help="Use this flag to omit any prompts", required=False)
    parser.add_argument("-O", "--optimize", action="store_true", help="Use to optimize DNA sequence instead translating it.", required=False)
    args = parser.parse_args()

    # Input files
    if args.table:
        codon_table = pandas.read_csv(args.table, header=None)
    else:
        codon_table = pandas.read_csv(path.join(path.dirname(path.abspath(__file__)), "data/sample_table.csv"), header=None)

    sequences = seqtools.io.fasta.open_fasta(args.input)

    # Protein translation
    if args.protein:
        solution = seqtools.seqtools.protein_to_dna(sequences, codon_table)
    elif not args.protein:
        solution = seqtools.seqtools.dna_operation(sequences, codon_table, args.force, args.optimize)

    # Saving/Printing solution(s)
    if args.output:
        path_save = path.join(path.join(getcwd(), args.output))
        msg_saved = "Output saved to `{0}`".format(path_save)
        seqtools.io.output.writer(solution, path_save)
        print(msg_saved)
    else:
        print()
        print(seqtools.io.output.writer(solution))

