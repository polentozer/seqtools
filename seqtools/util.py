import os
import sys
import pandas
from distutils.util import strtobool


def make_triplets(sequence):
    """
    Makes list of chunks 3 characters long from a sequence.
    """
    return [sequence[start:start + 3] for start in range(0, len(sequence), 3)]


def codon_table_parser(path_to_table):
    if not path_to_table:
        return None

    filenmae, extension = os.path.splitext(path_to_table)

    if extension == '.csv':
        return pandas.read_csv(path_to_table, header=None)

    raise ValueError('Codon table is not in the CSV format.')


def bool_user_prompt(question):
    '''
    Prompts user for YES/NO response
    '''
    sys.stdout.write(f'\n{question} [y/n]\n')
    while True:
        try:
            return strtobool(input().lower())
        except ValueError:
            sys.stdout.write('Please respond with "yes" or "no".\n')
