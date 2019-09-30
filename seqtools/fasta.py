import sys
from seqtools.modules import Nucleotide, Protein


def open_fasta(file_input_paths, protein=False):
    '''
    Function for opening `.fasta` files. Files can contain multiple sequences.
    Returns list of sequence objects.
    '''
    sequences = []
    seq_id = ''
    sequence = ''

    for input_path in file_input_paths:
        with open(input_path, 'r') as input_file:
            data = input_file.readlines()

        for line in data:
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


def write_fasta(solution_list, path=None):
    """
    Simple function for writing `.fasta` files.
    """
    if path:
        with open(path, 'w') as file_out:
            for sequence in solution_list:
                if sequence:
                    file_out.write(f'\n{sequence}\n')
    else:
        for sequence in solution_list:
            if sequence:
                sys.stdout.write(f'\n{sequence}\n')
