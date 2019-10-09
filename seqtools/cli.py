#cSpell: Disable#
'''
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
  seqtools trans -f input.fasta -o output.fasta
  seqtools trans -Of input.fasta -o output.fasta


Help:
  For help using this tool, please open an issue on the Github repository:
  https://github.com/polentozer/seqtools
'''
import os
import sys
import pandas
import logging
import argparse
import logging.config
from . import __version__ as VERSION
from seqtools.util import load_codon_table
from seqtools.generator import generate_dna
from seqtools.fasta import open_fasta, write_fasta
from seqtools.modules import Nucleotide, Protein, Restriction_Enzyme
from seqtools.seqtools_config import LOGGING_CONFIG, RESTRICTION_ENZYMES

### LOGGER ###
logging.config.dictConfig(LOGGING_CONFIG)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def cli_argparser():
    '''Argument parser for seqtools.'''

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
        'gen',
        help='Generate DNA sequence')

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
        '-e', '--restriction', action='store_true', default=False, required=False,
        help='Do not check for restriction sites')
    generate_parser.add_argument(
        '-r', '--ratio_gc', action='store_true', default=True, required=False,
        help='Do not limit the GC content between 0.4 and 0.6')

    # ==================
    # analyze sequence
    # ==================
    analyze_parser = subparser.add_parser(
        'anal',
        help='Analyze sequence')

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
        '-p', '--protein', action='store_true', required=False,
        help='Use this flag when input is a protein sequence')
    analyze_parser.add_argument(
        '-f', '--input_file', type=str, nargs='+', required=False,
        help='Path to input .fasta files')
    analyze_parser.add_argument(
        '-i', '--input', type=str, required=False,
        help='Paste raw sequence')

    # ==================
    # transform
    # ==================
    transform = subparser.add_parser(
        'trans',
        help='Transform/optimize DNA/protein sequence')

    # ------------------
    # transform options
    # ------------------
    table_input = transform.add_mutually_exclusive_group()
    transform.add_argument(
        '-O', '--optimize', action='store_true', default=False,
        help='Use this flag to change output to optimized DNA sequences')
    transform.add_argument(
        '-f', '--input_file', type=str, nargs='+', required=False,
        help='Path to input .fasta files')
    transform.add_argument(
        '-i', '--input', type=str, required=False,
        help='Paste raw sequence')
    transform.add_argument(
        '-o', '--output', type=str, nargs='?',
        help='Path for the output fasta file')
    transform.add_argument(
        '-P', '--protein', action='store_true', required=False, default=False,
        help='Use this flag when input is a protein sequence')
    transform.add_argument(
        '-M', '--maximum', action='store_true', required=False, default=False,
        help='Use this flag to use only the most common codons when optimizing sequences')
    transform.add_argument(
        '-t', '--type', type=str, required=False, default=None,
        help='Use this flag if you want your sequences transformed to parts of certain type')
    transform.add_argument(
        '-re', '--remove_cutsites', action='store_true', default=False,
        help='Use this flag when you want to remove cutsites from your DNA sequences. Used by default if "type" arguemnt is given')
    table_input.add_argument(
        '-on', '--organism_name', type=str, required=False, default=None,
        help='Organism name: ecoli/yeast/human/bsub/yali for spsum format')
    table_input.add_argument(
        '-oi', '--organism_id', type=int, required=False, default=None,
        help='Organism ID: taxonomy ID')

    # version and end of arguments
    parser.add_argument('--version', action='version', version=VERSION)

    return parser.parse_args()


def main():
    '''Main CLI entrypoint'''

    args = cli_argparser()

    if args.type or args.remove_cutsites:
        logger.debug('Generating RENZYME_LIST')
        renzyme_list = []
        for renzyme in RESTRICTION_ENZYMES.keys():
            renzyme_list.append(
                Restriction_Enzyme(
                    renzyme,
                    RESTRICTION_ENZYMES[renzyme]['substrate'],
                    RESTRICTION_ENZYMES[renzyme]['recognition'],
                    jump=RESTRICTION_ENZYMES[renzyme]['jump'],
                    overhang=RESTRICTION_ENZYMES[renzyme]['overhang'],
                    enzyme_description=RESTRICTION_ENZYMES[renzyme]['description']
                )
            )
    
    if args.commands == 'trans':
        logger.debug('SEQTOOLS-TRANSFORM')
        if args.organism_name or args.organism_id:
            logger.info(f'Loading codon table: species={args.organism_name} tax_id={args.organism_id}')
            codon_table = load_codon_table(species=args.organism_name, taxonomy_id=args.organism_id)
        else:
            logger.warning('Using default codon table: YALI')
            codon_table = load_codon_table(species='yali')
        if args.input_file:
            logger.info('Loading sequences...')
            sequences = open_fasta(args.input_file, protein=args.protein)
        else:
            if args.protein:
                sequences = [Protein('[SEQ_ID]', str(args.input))]
                logger.info('Reverse-translating proteins...')
                solution = [protein.reverse_translate(table=codon_table, maximum=args.maximum) for protein in sequences]
            else:
                sequences = [Nucleotide('[SEQ_ID]', str(args.input))]
                if args.optimize:
                    logger.info('Optimizing DNA sequences...')
                    solution = [dna.optimize_codon_usage(table=codon_table, maximum=args.maximum) for dna in sequences]
                else:
                    logger.info('Translating DNA sequences...')
                    solution = [dna.translate(table=codon_table) for dna in sequences]
                if args.type:
                    if args.type in ('1', '2', '3', '3a', '3b', '3t', '4', '4a', '4b', '5', '6', '7', '8', '8a', '8b'):
                        logger.info(f'Changing sequences to type{args.type} parts')
                        solution = [dna.remove_cutsites(renzyme_list, table=codon_table).make_part(part_type=args.type) for dna in solution]
                    else:
                        logger.warning(f'Part type{args.type} does not exist; skipping transformation into parts')
                elif args.remove_cutsites:
                    logger.info('Removing cutsites from DNA sequences...')
                    solution = [dna.remove_cutsites(renzyme_list, table=codon_table) for dna in solution]

        # saving/printing solutions
        if args.output:
            path_save = os.path.join(os.path.join(os.getcwd(), args.output))
            write_fasta(solution, path_save)
            logger.info(f'Output saved to `{path_save}`')
        else:
            for sequence in solution:
                sys.stdout.write(f'\n{sequence.fasta}\n')
        logger.info('DONE')

    elif args.commands == 'gen':
        logger.debug('SEQTOOLS-GENERATOR')
        generate_dna(
            args.length,
            args.homopolymer,
            args.gc_streach,
            args.restriction,
            args.ratio_gc
        )
        logger.info('DONE')

    elif args.commands == 'anal':
        logger.debug('SEQTOOLS-ANALYZER')
        sequences = open_fasta(args.input_file, protein=args.protein)
        logging.info('K-MER analysis')
        for seq in sequences:
            logging.info(f'{seq.kmer_occurrence(threshold=args.threshold, length=args.kmer)}')
        logger.info('DONE')

    else:
        logging.warning('No arguments were given. DONE')
