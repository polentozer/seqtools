###
##
#
from seqtools import __about__
import argparse


def cli_argparser():
    # arguments
    parser = argparse.ArgumentParser(description=__about__.__summary__)
    subparser = parser.add_subparsers(dest='commands')

    # generate dna
    generate_parser = subparser.add_parser("gen", help="Generate DNA sequence")

    # generator options
    generate_parser.add_argument('length', action='store', type=int, help='How long should generated DNA be')
    generate_parser.add_argument('-s', '--short', action='store', default=4, required=False, type=int, help='Number of allowed single letter repeats')
    generate_parser.add_argument('-l', '--long', action='store', default=6, required=False, type=int, help=' Maximum allowed GC/AT stretch')
    generate_parser.add_argument('-t', '--type2', action='store_true', default=False, required=False, help='Do not check for restriction sites')
    generate_parser.add_argument('-g', '--gc', action='store_true', default=False, required=False, help='Do not limit the GC content between 0.4 and 0.6')

    # translate
    translate_parser = subparser.add_parser("trans", help="Translate/optimize DNA/protein sequence")

    # translate options
    group = translate_parser.add_mutually_exclusive_group()
    group.add_argument("-A", "--analyze", action="store_true", help="Use this flag to perform analysis on your sequences")
    group.add_argument("-O", "--optimize", action="store_true", help="Use this flag to optimize DNA sequence instead translating it.")
    translate_parser.add_argument("-f", "--force", action="store_true", help="Use this flag to omit any prompts and force the process through", required=False)
    translate_parser.add_argument("-i", "--input", help="Path to input 'fasta' files", type=str, nargs='+', required=True)
    translate_parser.add_argument("-o", "--output", help="Path for the output fasta file", type=str, nargs='?')
    translate_parser.add_argument("-p", "--protein", action="store_true", help="Use this flag when working with protein sequences", required=False)
    translate_parser.add_argument("-t", "--table", help="Path to codon usage table in csv format: 'aminoacid,triplet,value')", type=str, required=False)

    # version and end of arguments
    parser.add_argument("-V", "--version", action="version", version=__about__.__version__)

    return parser.parse_args()
