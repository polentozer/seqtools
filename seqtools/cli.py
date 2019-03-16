"""
seqtools

Usage:
  seqtools gen 20
  seqtools trans -i input.fasta
  seqtools -h | --help
  seqtools --version

Options:
  -h --help                         Show this screen.
  --version                         Show version.

Examples:
  seqtools gen 20
  seqtools gen 45 -
  seqtools trans -i input.fasta
  seqtools trans -i input.fasta -o output.fasta
  seqtools trans -O -i input.fasta -o output.fasta


Help:
  For help using this tool, please open an issue on the Github repository:
  https://github.com/polentozer/seqtools
"""

import os
import sys
import pandas
import argparse
from . import __version__ as VERSION
from seqtools.util import codon_table_parser
from seqtools.io.fasta import open_fasta, write_fasta
from seqtools.generate.generator import generate_dna


def cli_argparser():
    """Argument parser for seqtools."""

    # arguments
    parser = argparse.ArgumentParser(description='A command line tool for manipulating `.fasta` files')
    subparser = parser.add_subparsers(dest='commands')

    # generate dna
    generate_parser = subparser.add_parser("gen", help="Generate DNA sequence")

    # generator options
    generate_parser.add_argument('length', action='store', type=int, help='How long should generated DNA be')
    generate_parser.add_argument('-n', '--single_repeats', action='store', default=4, required=False, type=int, help='Number of allowed single nucleotide repeats')
    generate_parser.add_argument('-g', '--gc_streach', action='store', default=6, required=False, type=int, help=' Maximum allowed GC/AT stretch')
    generate_parser.add_argument('-t', '--type2', action='store_true', default=False, required=False, help='Do not check for restriction sites')
    generate_parser.add_argument('-r', '--ratio_gc', action='store_true', default=False, required=False, help='Do not limit the GC content between 0.4 and 0.6')

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
    parser.add_argument("--version", action="version", version=VERSION)

    return parser.parse_args()


def main():
    """Main CLI entrypoint"""

    args = cli_argparser()

    # Old simple stuff
    if args.commands == 'trans':

        # Input DATA
        if args.table:
            codon_table = codon_table_parser(args.table)
        else:
            sys.stdout.write("\n### Using sample codon usage table!!! ###\n")
            codon_table = pandas.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/sample_table.csv"), header=None)

        # protein translation
        if args.protein and not args.analyze:
            sequences = open_fasta(args.input, protein=True)
            solution = [protein.reverse_translate(codon_table) for protein in sequences]
        else:
            sequences = open_fasta(args.input)
            if args.optimize:
                solution = [dna.optimize_codon_usage(codon_table) for dna in sequences]
            elif args.analyze:
                raise NotImplemented()
            else:
                solution = [dna.translate(codon_table) for dna in sequences]

        # saving/printing solutions
        if args.output:
            path_save = os.path.join(os.path.join(os.getcwd(), args.output))
            msg_saved = "Output saved to `{0}`".format(path_save)
            write_fasta(solution, path_save)
            sys.stdout.write(msg_saved)
        else:
            write_fasta(solution)

    elif args.commands == 'gen':
        generate_dna(args.length, args.single_repeats, args.gc_streach, args.type2, args.ratio_gc)

    else:
        sys.stdout.write('No arguments were given.')
