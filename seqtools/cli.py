#cSpell: Disable#
"""
seqtools

Usage:
  seqtools gen 20
  seqtools trans -i input.fasta
  seqtools -h | --help
  seqtools --version

Options:
  -h --help                         Show this screen
  --version                         Show version
  gen                               Generate DNA sequence
  trans                             Transform sequence(s)
  anal                              Analyze sequences

Examples:
  seqtools gen 20
  seqtools gen 45 -n 5 -g 10
  seqtools trans -f input.fasta
  seqtools trans -Oi ATGATATG
  seqtools trans -Vf input.fasta
  seqtools trans -f input.fasta -o output.fasta
  seqtools trans -Of input.fasta -o output.fasta


Help:
  For help using this tool, please open an issue on the Github repository:
  https://github.com/polentozer/seqtools
"""

import os
import sys
import pandas
import argparse
from . import __version__ as VERSION
from seqtools.generator import generate_dna
from seqtools.util import codon_table_parser
from seqtools.modules import Nucleotide, Protein
from seqtools.fasta import open_fasta, write_fasta


def cli_argparser():
    """Argument parser for seqtools."""

    ## ===========
    ## arguments
    ## ===========
    parser = argparse.ArgumentParser(
        description='A command line tool for manipulating `.fasta` files')

    subparser = parser.add_subparsers(dest='commands')

    # ==================
    # generate dna
    # ==================
    generate_parser = subparser.add_parser(
        "gen",
        help="Generate DNA sequence")

    # ------------------
    # generator options
    # ------------------
    generate_parser.add_argument(
        'length', action='store', type=int,
        help='How long should generated DNA be')
    generate_parser.add_argument(
        '-n', '--homopolymer', action='store', default=4, required=False, type=int,
        help='Number of allowed single nucleotide repeats')
    generate_parser.add_argument(
        '-g', '--gc_streach', action='store', default=6, required=False, type=int,
        help=' Maximum allowed GC/AT stretch')
    generate_parser.add_argument(
        '-t', '--type2', action='store_true', default=False, required=False,
        help='Do not check for restriction sites')
    generate_parser.add_argument(
        '-r', '--ratio_gc', action='store_true', default=False, required=False,
        help='Do not limit the GC content between 0.4 and 0.6')

    # ==================
    # analyze sequence
    # ==================
    analyze_parser = subparser.add_parser(
        "anal",
        help="Analyze sequence")

    # ------------------
    # analyze options
    # ------------------
    analyze_parser.add_argument(
        '-t', '--threshold', action='store', type=int,
        help='How many times should a kmer be repeated to display it')
    analyze_parser.add_argument(
        '-k', '--kmer', action='store', type=int,
        help='Length of kmers')
    analyze_parser.add_argument(
        "-p", "--protein", action="store_true", required=False,
        help="Use this flag when input is a protein sequence")
    analyze_parser.add_argument(
        "-f", "--input_file", type=str, nargs='+', required=False,
        help="Path to input 'fasta' files")
    analyze_parser.add_argument(
        "-i", "--input", type=str, required=False,
        help="Paste raw sequence")

    # ==================
    # translate
    # ==================
    translate_parser = subparser.add_parser(
        "trans",
        help="Translate/optimize DNA/protein sequence")

    # ------------------
    # translate options
    # ------------------
    group = translate_parser.add_mutually_exclusive_group()
    group.add_argument(
        "-V", "-optvalue", action="store_true", dest="optvalue",
        help="Use this flag to perform optimization value analysis on your sequences")
    group.add_argument(
        "-O", "--optimize", action="store_true",
        help="Use this flag to optimize DNA sequence instead translating it.")
    translate_parser.add_argument(
        "-f", "--input_file", type=str, nargs='+', required=False,
        help="Path to input 'fasta' files")
    translate_parser.add_argument(
        "-i", "--input", type=str, required=False,
        help="Paste raw sequence")
    translate_parser.add_argument(
        "-o", "--output", type=str, nargs='?',
        help="Path for the output fasta file")
    translate_parser.add_argument(
        "-p", "--protein", action="store_true", required=False,
        help="Use this flag when input is a protein sequence")
    translate_parser.add_argument(
        "-t", "--table", type=str, required=False,
        help="Path to codon usage table in csv format: 'aminoacid,triplet,value'")

    # ------------------
    # version and end of arguments
    # ------------------
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
            print("\n### Using sample codon usage table!!! ###\n")
            codon_table = pandas.read_csv(os.path.join(
                os.path.dirname(os.path.abspath(__file__)), "data/sample_table.csv"),
                header=None)

        # protein translation
        if args.protein:
            if args.input_file:
                sequences = open_fasta(args.input_file, protein=True)
            else:
                sequences = [Protein('0x0', str(args.input))]
            solution = [protein.reverse_translate(codon_table) for protein in sequences]

        else:
            if args.input_file:
                sequences = open_fasta(args.input_file)
            else:
                sequences = [Nucleotide('0x0', str(args.input))]

            if args.optimize:
                solution = [dna.optimize_codon_usage(codon_table) for dna in sequences]
            elif args.optvalue:
                solution = []
                print('\nOptimization values are:')
                for dna in sequences:
                    sys.stdout.write(
                        f'\n{dna.sequence_id}: {dna.optimization_value(codon_table)}\n')
            else:
                solution = [dna.translate(codon_table) for dna in sequences]

        # saving/printing solutions
        if args.output:
            path_save = os.path.join(os.path.join(os.getcwd(), args.output))
            msg_saved = "Output saved to `{0}`".format(path_save)
            write_fasta(solution, path_save)
            print(msg_saved)
        else:
            write_fasta(solution)

    elif args.commands == 'gen':
        generate_dna(
            args.length,
            args.homopolymer,
            args.gc_streach,
            args.type2,
            args.ratio_gc
        )

    elif args.commands == 'anal':
        sequences = open_fasta(args.input_file, protein=args.protein)
        print('\nK-MER analysis\n')
        for seq in sequences:
            print(f'{seq.kmer_occurrence(threshold=args.threshold, length=args.kmer)}\n')

    else:
        print('No arguments were given.')
