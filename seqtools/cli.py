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
import argparse
import logging
import logging.config
from . import __version__ as VERSION
from seqtools.seqtools_config import LOGGING_CONFIG, RESTRICTION_ENZYMES

### LOGGER ###
# TODO: Make argumental logging levels -v (verbose) -w (debug)
logging.config.dictConfig(LOGGING_CONFIG)
logger = logging.getLogger(__name__)

from seqtools.util import load_codon_table
from seqtools.generator import generate_dna
from seqtools.fasta import open_fasta, write_fasta
from seqtools.modules import Nucleotide, Protein, Restriction_Enzyme


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
        '-n', '--homopolymer', action='store', default=10, required=False, type=int,
        help='Number of allowed single nucleotide repeats')
    generate_parser.add_argument(
        '-g', '--gc_streatch', action='store', default=20, required=False, type=int,
        help=' Maximum allowed GC/AT stretch')
    generate_parser.add_argument(
        '-e', '--restriction', action='store_true', default=False, required=False,
        help='Do not check for restriction sites')
    generate_parser.add_argument(
        '-r', '--ratio_gc', action='store_true', default=True, required=False,
        help='Do not limit the GC content between 0.4 and 0.6')
    generate_parser.add_argument(
        '-P', '--protein', action='store_true', default=False, required=False,
        help='Use this flag to generate a random protein sequence of chosen length')

    # ==================
    # analyze sequence
    # ==================
    analyze_parser = subparser.add_parser(
        'anal',
        help='Analyze sequence')

    # ------------------
    # analyze options
    # ------------------
    analyze_sequence_input = analyze_parser.add_mutually_exclusive_group()
    analyze_parser.add_argument(
        '-t', '--threshold', action='store', type=int, required=False,
        help='How many times should a kmer be repeated to display it')
    analyze_parser.add_argument(
        '-k', '--kmer', action='store', type=int, required=False,
        help='Length of kmers')
    analyze_parser.add_argument(
        '-p', '--protein', action='store_true', required=False,
        help='Use this flag when input is a protein sequence')
    analyze_sequence_input.add_argument(
        '-f', '--input_file', type=str, nargs='+', required=False,
        help='Path to input .fasta files')
    analyze_sequence_input.add_argument(
        '-i', '--input', type=str, required=False,
        help='Paste raw sequence')
    analyze_parser.add_argument(
        '-g', '--graph', action='store_true', required=False, default=False,
        help='Use this flag to graph codon usage frequency')


    # ==================
    # transform
    # ==================
    transform = subparser.add_parser(
        'trans',
        help='Transform/optimize DNA/protein sequence')

    # ------------------
    # transform options
    # ------------------
    sequence_input = transform.add_mutually_exclusive_group()
    table_input = transform.add_mutually_exclusive_group()
    maniputaion = transform.add_mutually_exclusive_group()
    graph = transform.add_mutually_exclusive_group()
    maniputaion.add_argument(
        '-O', '--optimize', action='store_true', default=False,
        help='Use this flag to change output to optimized DNA sequences')
    maniputaion.add_argument( # TODO
        '-SO', '--sourceoptimize', type=int, default=None,
        help='Use this flag to optimize sequence based on codons usage from source organism (tax id)')
    maniputaion.add_argument(
        '-T', '--translate', action='store_true', default=False,
        help='Use this flag to change output to proteins translated from given DNAs')
    sequence_input.add_argument(
        '-f', '--input_file', type=str, nargs='+', required=False,
        help='Path to input .fasta files')
    sequence_input.add_argument(
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
    graph.add_argument(
        '-g', '--graph', action='store_true', required=False, default=False,
        help='Use this flag to show codon usage graph')
    graph.add_argument(
        '-G', '--graph_optimized', action='store_true', required=False, default=False,
        help='Use this flag to show codon usage graph for optimized and basic CDS')


    # ==================
    # version and end of arguments
    # ==================
    parser.add_argument('--version', action='version', version=VERSION)

    return parser.parse_args()


def main():
    '''Main CLI entrypoint'''
    args = cli_argparser()
    
    # TRANSFORMATIONS
    if args.commands == 'trans':
        logger.debug('SEQTOOLS-TRANSFORM')

        sequences, solution = [], []

        # Codon table input
        if args.organism_name or args.organism_id:
            logger.info(f'Loading codon table: species={args.organism_name} tax_id={args.organism_id}')
            codon_table = load_codon_table(species=args.organism_name, taxonomy_id=args.organism_id)
        else:
            logger.warning('Using default codon table: YALI')
            codon_table = load_codon_table(species='yali')

        # Sequences input
        if args.input_file:
            logger.info('Loading sequences...')
            sequences = open_fasta(args.input_file, protein=args.protein)
        else:
            if args.protein:
                sequences = [Protein('[SEQ_ID]', str(args.input))]
                logger.info('Reverse-translating proteins...')
            else:
                sequences = [Nucleotide('[SEQ_ID]', str(args.input))]

        # Protein reverse translation
        if args.protein and sequences:
            solution = [protein.reverse_translate(table=codon_table, maximum=args.maximum) for protein in sequences]

        # DNA sequences
        elif sequences:
            # Optimization
            if args.optimize:
                logger.info('Optimizing DNA sequences...')
                solution = [dna.optimize_codon_usage(table=codon_table, maximum=args.maximum) for dna in sequences]
            # Translation
            elif args.translate:
                logger.info('Translating DNA sequences...')
                solution = [dna.translate(table=codon_table) for dna in sequences]
            # Special optimization TODO:
            elif args.sourceoptimize:
                source_organism_table = load_codon_table(taxonomy_id=args.sourceoptimize)
                solution = [dna.source_optimize(source_table=source_organism_table, table=codon_table) for dna in sequences]
            
            # Restriction enzyme recognition sites removal
            if args.remove_cutsites:
                logger.debug('Generating RENZYME_LIST')

                # Restriction enzyme parser for RESTRICTION_ENZYME dictionary (found in seqtools_config.py)
                renzyme_list = []
                for renzyme in RESTRICTION_ENZYMES.keys():
                    renzyme_list.append(
                        Restriction_Enzyme(
                            renzyme,
                            RESTRICTION_ENZYMES[renzyme]['substrate'],
                            RESTRICTION_ENZYMES[renzyme]['recognition'],
                            jump=RESTRICTION_ENZYMES[renzyme]['jump'],
                            overhang=RESTRICTION_ENZYMES[renzyme]['overhang'],
                            enzyme_description=RESTRICTION_ENZYMES[renzyme]['description']))

                if solution:
                    target = solution
                else:
                    target = sequences
                logger.info('Removing cutsites from DNA sequences...')
                solution = [dna.remove_cutsites(renzyme_list, table=codon_table) for dna in target]

            # Appending prefixes and suffixes for GGE parts
            if args.type:
                if args.type in ('1', '2', '3', '3a', '3b', '3t', '4', '4a', '4b', '5', '6', '7', '8', '8a', '8b'):
                    if args.graph or args.graph_optimized:
                        storage = solution
                    if solution:
                        target = solution
                    else:
                        target = sequences
                    logger.info(f'Changing sequences to type{args.type} parts')
                    solution = [dna.make_part(part_type=args.type) for dna in target]
                else:
                    logger.warning(f'Part type{args.type} does not exist; skipping transformation into parts')

        # No sequences
        else:
            logger.error('Fasta file is empty or no sequence provided...')

        # Saving/printing solutions
        if args.output:
            path_save = os.path.join(os.path.join(os.getcwd(), args.output))
            write_fasta(solution, path_save)
            logger.info(f'Output saved to `{path_save}`')
        else:
            for sequence in solution:
                sys.stdout.write(f'\n{sequence.fasta}\n')
        
        # Drawing graphs
        if args.graph:
            logger.info('Drawing graphs...')
            if not args.type:
                storage = solution
            for seq in storage:
                seq.graph_codon_usage(window=10, table=codon_table)
        elif args.graph_optimized:
            logger.info('Drawing graphs...')
            for seq_original, seq_optimized in zip(sequences, storage):
                seq_optimized.graph_codon_usage(window=10, other=seq_original, table=codon_table)
        
        logger.info('DONE')

    # GENERATOR
    elif args.commands == 'gen':
        logger.debug('SEQTOOLS-GENERATOR')
        generate_dna(
            args.length,
            homopolymer=args.homopolymer,
            gc_stretch=args.gc_streatch,
            restriction=args.restriction,
            ratio_gc=args.ratio_gc,
            protein=args.protein
        )
        logger.info('DONE')

    # ANALYZER
    elif args.commands == 'anal':
        logger.debug('SEQTOOLS-ANALYZER')
        sequences = open_fasta(args.input_file, protein=args.protein)
        logging.info('K-MER analysis')
        if sequences:
            for seq in sequences:
                logging.info(f'{seq.kmer_analysis(threshold=args.threshold, length=args.kmer)}')
        else:
            logger.error('Fasta file is empty or no sequence provided...')
        logger.info('DONE')

    # NO ARGUMENTS
    else:
        logging.warning('No arguments given.\nDONE')
