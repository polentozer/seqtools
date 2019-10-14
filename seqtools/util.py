#cSpell: Disable#
import os
import signal
import pandas
from random import random
from seqtools.cli import logger
from distutils.util import strtobool
from seqtools.seqtools_config import LOGGING_CONFIG


CODON_USAGE_DB = f'{os.path.dirname(__file__)}/data/codon_usage.spsum'
CUSTOM_CODON_USAGE_DB = f'{os.path.dirname(__file__)}/data/custom_table.spsum'

COMMON_SPECIES = {
    'ecoli': '83333',
    'yeast':  '4932',
    'human': '9606',
    'bsub': '1432',
    'yali': '284591'}

CODONS = [
    'CGA', 'CGC', 'CGG', 'CGT', 'AGA', 'AGG', 'CTA', 'CTC',
    'CTG', 'CTT', 'TTA', 'TTG', 'TCA', 'TCC', 'TCG', 'TCT',
    'AGC', 'AGT', 'ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'CCC',
    'CCG', 'CCT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC',
    'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'AAA', 'AAG',
    'AAC', 'AAT', 'CAA', 'CAG', 'CAC', 'CAT', 'GAA', 'GAG',
    'GAC', 'GAT', 'TAC', 'TAT', 'TGC', 'TGT', 'TTC', 'TTT',
    'ATA', 'ATC', 'ATT', 'ATG', 'TGG', 'TAA', 'TAG', 'TGA']

STANDARD_GENETIC_CODE = [
    'R', 'R', 'R', 'R', 'R', 'R', 'L', 'L',
    'L', 'L', 'L', 'L', 'S', 'S', 'S', 'S',
    'S', 'S', 'T', 'T', 'T', 'T', 'P', 'P',
    'P', 'P', 'A', 'A', 'A', 'A', 'G', 'G',
    'G', 'G', 'V', 'V', 'V', 'V', 'K', 'K',
    'N', 'N', 'Q', 'Q', 'H', 'H', 'E', 'E',
    'D', 'D', 'Y', 'Y', 'C', 'C', 'F', 'F',
    'I', 'I', 'I', 'M', 'W', '*', '*', '*']

####### TIMEOUT DECORATOR #######
class TimeoutError(Exception):
    def __init__(self, value='Timed Out'):
        self.value = value
    def __str__(self):
        return repr(self.value)

def timeout(seconds_before_timeout):
    def decorate(f):
        def handler(signum, frame):
            raise TimeoutError()
        def new_f(*args, **kwargs):
            old = signal.signal(signal.SIGALRM, handler)
            signal.alarm(seconds_before_timeout)
            try:
                result = f(*args, **kwargs)
            finally:
                # reinstall the old signal handler
                signal.signal(signal.SIGALRM, old)
                # cancel the alarm
                signal.alarm(0)
            return result
        new_f.__name__ = f.__name__
        return new_f
    return decorate
#################################

def bool_user_prompt(nucleotide_obj, operation):
    '''Prompts user for YES/NO response. Version 2.'''
    logger.error(f'Sequence with ID {nucleotide_obj.sequence_id} is not a CDS. {operation} anyway? [y/n]')
    while True:
        try:
            if strtobool(input().lower()):
                logger.debug(f'Responded with TRUE')
                nucleotide_obj.sequence_id += '|FORCED'
            return nucleotide_obj
        except ValueError:
            logger.warning('Please respond with "yes" or "no".')


def sequence_match(string, search):
    '''Returns TRUE if sequence matches condition in search'''
    logger.debug('Sequence matching...')
    return not bool(search(string))


def load_codon_table(species=None, taxonomy_id=None, custom=False):
    '''Load a codon table based on the organism's species ID'''
    logger.debug(f'load_codon_table(species={species}, taxonomy_id={taxonomy_id}, custom={custom})')

    if species in COMMON_SPECIES:
        taxonomy_id = COMMON_SPECIES[species]

    if custom:
        codon_usage = CUSTOM_CODON_USAGE_DB
    else:
        codon_usage = CODON_USAGE_DB

    with open(codon_usage) as cu:
        for header in cu:
            codon_counts = cu.readline()

            taxid, species = header.strip().split(':')[:2]

            if taxonomy_id:
                taxonomy_id = str(taxonomy_id)

            if taxonomy_id and taxonomy_id != taxid:
                continue

            table = list(
                zip(CODONS, STANDARD_GENETIC_CODE, [int(x) for x in codon_counts.split()]))
            table = pandas.DataFrame(table, columns=['Triplet', 'AA', 'Number'])
            table.set_index(['AA', 'Triplet'], inplace=True)
            table.sort_index(inplace=True)

            table['Fraction'] = table.groupby('AA').transform(lambda x: x / x.sum())
            break

    return table


def get_codon(codons, maximum=False, recode=False, skip=[]):
    '''Returns a "locally-optimized" codon. Locally-optimized = mimics the
    codon frequency in the table. Maximum uses the most common codon.'''
    logger.debug(f'get_codon({list(codons.index)}, maximum={maximum}, recode={recode}, skip={skip})')
    if recode:
        codons = codons.loc[[cod for cod in codons.index if cod not in skip]]
    # Deterministic allocation of codon based on the highest frequency
    if maximum:
        return codons.Fraction.idxmax()
    # Stochastic allocation of codon
    x = codons.Fraction.cumsum() / codons.Fraction.cumsum().max() < random()

    return codons.iloc[x.sum()].name


def codon_table_10plus(table):
    '''Return a codon table only representing codons with > 10% occurrence frequency.'''
    logger.debug('Generating codon table 10+...')
    table = table[table.Fraction >= 0.1]
    table = table.groupby(level=0).transform(lambda x: x / x.sum())

    return table
