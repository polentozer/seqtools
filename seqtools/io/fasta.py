from seqtools.modules import Nucleotide, Protein


def open_fasta(file_input_paths, protein=False):
    '''
    Function for opening `.fasta` files. Files can contain multiple sequences.
    Returns list of sequences.
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

    return sequences
