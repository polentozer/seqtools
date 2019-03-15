import os
import sys
import pandas
from seqtools.io.output import writer
from seqtools.cli import cli_argparser
from seqtools.io.fasta import open_fasta
from seqtools.util import codon_table_parser
from seqtools.generate.generator import generate_dna


def main():

    args = cli_argparser()

    # Old simple stuff
    if args.commands == 'trans':

        # Input DATA
        if args.table:
            codon_table = codon_table_parser(args.table)
        else:
            sys.stdout.write("\n### Using sample codon usage table!!! ###")
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
            # seqtools.dna_operation(sequences, codon_table, args.force, args.optimize, args.analyze)

        # saving/printing solutions
        if args.output:
            path_save = os.path.join(os.path.join(getcwd(), args.output))
            msg_saved = "Output saved to `{0}`".format(path_save)
            writer(solution, path_save)
            sys.stdout.write(msg_saved)
        else:
            writer(solution)

    elif args.commands == 'gen':
        generate_dna(args.length, args.short, args.long, args.type2, args.gc)

    else:
        sys.stdout.write('No arguments were given.')
