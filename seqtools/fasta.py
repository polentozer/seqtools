import logging
from seqtools.modules import Nucleotide, Protein


def open_fasta(file_input_paths, protein=False):
    '''Function for opening `.fasta` files. Files can contain multiple sequences.
    Returns list of sequence objects.'''

    logger = logging.getLogger(__name__)

    logger.info('Opening fasta file(s)...')
    logger.debug(f'open_fasta({file_input_paths}, protein={protein})')
    sequences = []
    seq_id, sequence = '', ''

    for input_path in file_input_paths:
        with open(input_path, 'r') as input_file:
            for line in input_file.readlines():
                temp = line.strip()
                if line[0] == '>':
                    if seq_id != '':
                        if protein:
                            sequences.append(Protein(seq_id, sequence))
                        else:
                            sequences.append(Nucleotide(seq_id, sequence))
                    seq_id = temp[:40]
                    sequence = ''
                else:
                    sequence += temp.upper()
            else:
                if protein:
                    sequences.append(Protein(seq_id, sequence))
                else:
                    sequences.append(Nucleotide(seq_id, sequence))

    return sequences


def write_fasta(sequence_list, path=''):
    '''Simple function for writing `.fasta` files.'''
    with open(f'{path}output.fasta', 'w') as file_out:
        for sequence in sequence_list:
            if sequence:
                file_out.write(sequence.fasta)
