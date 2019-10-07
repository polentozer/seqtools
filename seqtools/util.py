#cSpell: Disable#
import os
import sys
import pandas
import numpy as np
from distutils.util import strtobool


def codon_table_parser(path_to_table):
    if not path_to_table:
        return None

    _, extension = os.path.splitext(path_to_table)

    if extension == '.csv':
        return pandas.read_csv(path_to_table, header=None)

    raise ValueError('Codon table is not in the CSV format.')


def bool_user_prompt(question):
    '''Prompts user for YES/NO response'''
    print(f'\n{question} [y/n]')
    while True:
        try:
            return strtobool(input().lower())
        except ValueError:
            print('Please respond with "yes" or "no".')


def sequence_match(string, search):
    '''Returns TRUE if sequence matches condition in search'''
    return not bool(search(string))


def load_codon_table(species=None, taxonomy_id=None, custom=False):
    '''Load a codon table based on the organism's species ID'''
    CODON_USAGE_DB = f'{os.path.dirname(__file__)}/data/codon_usage.spsum'
    CUSTOM_CODON_USAGE_DB = f'{os.path.dirname(__file__)}/data/custom_table.spsum'
    COMMON_SPECIES = {
        'ecoli': '83333',
        'yeast':  '4932',
        'human': '9606',
        'bsub': '1432',
        'yali': '284591'
    }

    codons = [
        'CGA', 'CGC', 'CGG', 'CGT',
        'AGA', 'AGG', 'CTA', 'CTC',
        'CTG', 'CTT', 'TTA', 'TTG',
        'TCA', 'TCC', 'TCG', 'TCT',
        'AGC', 'AGT', 'ACA', 'ACC',
        'ACG', 'ACT', 'CCA', 'CCC',
        'CCG', 'CCT', 'GCA', 'GCC',
        'GCG', 'GCT', 'GGA', 'GGC',
        'GGG', 'GGT', 'GTA', 'GTC',
        'GTG', 'GTT', 'AAA', 'AAG',
        'AAC', 'AAT', 'CAA', 'CAG',
        'CAC', 'CAT', 'GAA', 'GAG',
        'GAC', 'GAT', 'TAC', 'TAT',
        'TGC', 'TGT', 'TTC', 'TTT',
        'ATA', 'ATC', 'ATT', 'ATG',
        'TGG', 'TAA', 'TAG', 'TGA'
    ]
    standard_genetic_code = [
        'R', 'R', 'R', 'R',
        'R', 'R', 'L', 'L',
        'L', 'L', 'L', 'L',
        'S', 'S', 'S', 'S',
        'S', 'S', 'T', 'T',
        'T', 'T', 'P', 'P',
        'P', 'P', 'A', 'A',
        'A', 'A', 'G', 'G',
        'G', 'G', 'V', 'V',
        'V', 'V', 'K', 'K',
        'N', 'N', 'Q', 'Q',
        'H', 'H', 'E', 'E',
        'D', 'D', 'Y', 'Y',
        'C', 'C', 'F', 'F',
        'I', 'I', 'I', 'M',
        'W', '*', '*', '*']

    if species in COMMON_SPECIES:
        taxonomy_id = COMMON_SPECIES[species]

    if custom:
        codon_usage = CUSTOM_CODON_USAGE_DB
    else:
        codon_usage = CODON_USAGE_DB

    with open(codon_usage) as cu:
        for header in cu:
            codon_counts = cu.readline()

            taxid, species, _ = header.strip().split(':')[:3]

            if taxonomy_id and taxonomy_id != taxid:
                continue

            table = list(zip(codons, standard_genetic_code, [int(x) for x in codon_counts.split()]))
            table = pandas.DataFrame(table, columns=['Triplet', 'AA', 'Number'])
            table.set_index(['AA', 'Triplet'], inplace=True)
            table.sort_index(inplace=True)

            table['Fraction'] = table.groupby('AA').transform(lambda x: x / x.sum())
            break
    
    return table


def get_codon(table, amino, maximum=False):
    '''Returns a "locally-optimized" codon for a given amino acid
    based on the single letter code. Locally-optimized = mimics the
    codon frequency in the table. Maximum uses the most common codon.'''
    c = table.loc[amino]
    if maximum:
        return c[c.Fraction == c.Fraction.max()].index[0]
    return c.iloc[(c.Fraction.cumsum() / c.Fraction.cumsum().max() < np.random.rand()).sum()].name
