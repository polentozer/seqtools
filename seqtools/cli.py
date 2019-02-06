###
##
#

import argparse
import pandas
import seqtools
from os import path, getcwd


def main():
    """
    Main function with argument parsing.
    """

    # Arguments
    parser = argparse.ArgumentParser(description=seqtools.__summary__)
    subparser = parser.add_subparsers(dest='commands')

    ## generate dna
    generate_parser = subparser.add_parser("gen", help="Generate DNA sequence")

    ### generator options
    generate_parser.add_argument('length', action='store', type=int, help='How long should generated DNA be')
    generate_parser.add_argument('-s', '--short', action='store', default=4, required=False, type=int, help='Number of allowed single letter repeats')
    generate_parser.add_argument('-l', '--long', action='store', default=6, required=False, type=int, help=' Maximum allowed GC/AT stretch')
    generate_parser.add_argument('-t', '--type2', action='store_true', default=False, required=False, help='Do not check for restriction sites')
    generate_parser.add_argument('-g', '--gc', action='store_true', default=False, required=False, help='Do not limit the GC content between 0.4 and 0.6')

    ## translate
    translate_parser = subparser.add_parser("trans", help="Translate/optimize DNA/protein sequence")

    ### translate options
    group = translate_parser.add_mutually_exclusive_group()
    group.add_argument("-A", "--analyze", action="store_true", help="Use this flag to perform analysis on your sequences")
    group.add_argument("-O", "--optimize", action="store_true", help="Use this flag to optimize DNA sequence instead translating it.")
    translate_parser.add_argument("-f", "--force", action="store_true", help="Use this flag to omit any prompts and force the process through", required=False)
    translate_parser.add_argument("-i", "--input", help="Path to input 'fasta' files", type=str, nargs='+', required=True)
    translate_parser.add_argument("-o",  "--output", help="Path for the output fasta file", type=str, nargs='?')
    translate_parser.add_argument("-p", "--protein", action="store_true", help="Use this flag when working with protein sequences", required=False)
    translate_parser.add_argument("-t", "--table",help="Path to codon usage table in csv format: 'aminoacid,triplet,value')", type=str, required=False)

    ## Version and end of arguments
    parser.add_argument("-V", "--version", action="version", version=seqtools.__version__)
    args = parser.parse_args()


    # Old simple stuff
    if args.commands == 'trans':

        ## Input DATA
        if args.table:
            codon_table = pandas.read_csv(args.table, header=None)
        else:
            print("\n### Using sample codon usage table!!! ###")
            codon_table = pandas.read_csv(path.join(path.dirname(path.abspath(__file__)), "data/sample_table.csv"), header=None)

        ## Input FASTA file
        sequences = seqtools.io.fasta.open_fasta(args.input)


        ## Protein translation
        if args.protein and not args.analyze:
            solution = seqtools.seqtools.protein_to_dna(sequences, codon_table)
        else:
            solution = seqtools.seqtools.dna_operation(sequences, codon_table, args.force, args.optimize, args.analyze)

        ## Saving/Printing solution(s)
        if args.output:
            path_save = path.join(path.join(getcwd(), args.output))
            msg_saved = "Output saved to `{0}`".format(path_save)
            seqtools.io.output.writer(solution, path_save)
            print(msg_saved)
        else:
            print()
            print(seqtools.io.output.writer(solution))
    
    elif args.commands == 'gen':
        seqtools.generate.generate_dna.generate_dna(args.length, args.short, args.long, args.type2, args.gc)


# Python script suffix
if __name__ == "__main__":
    main()
