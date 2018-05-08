import pandas


def make_triplets(string):
    """
    Makes chunks of 3 characters from a long string.
    """
    return [string[start:start+3] for start in range(0, len(string), 3)]


def prompt_force(string):
    cont_prompt = input("Sequence with ID `{0}` is not a CDS, optimize anyway? (Y/n): ".format(string))
    
    if cont_prompt.upper() == "Y" or cont_prompt.upper() == "YES" or cont_prompt == "":
        return True
    else:
        return False


def dna_operation(sequences, codon_table, force, optimize=None):
    """
    Optimizes DNA triplets to use most common codons from the codon_table.
    """
    solutions = {}

    if optimize:
        for sequence_id, sequence in sequences.items():
            solutions[sequence_id] = codon_optimize(sequence_id, sequence, codon_table, force)
    else:
        for sequence_id, sequence in sequences.items():
            solutions[sequence_id] = translate_dna(sequence_id, sequence, codon_table, force)
    
    return solutions


def protein_to_dna(sequences, codon_table):
    """
    Translates protein to dna with a call to reverse_translate() function.
    """
    reverse_translations = {}

    for sequence_id, sequence in sequences.items():
        reverse_translations[sequence_id] = ["", None]

        for amino in sequence:
            reverse_translations[sequence_id][0] += reverse_translate(amino, codon_table)
    
    return reverse_translations


def reverse_translate(amino, codon_table):
    """
    Reverse-translates aminoacid back to the most used codon from codon table.
    """
    temp_values = []
    temp_triplets = []

    for i, row in codon_table.loc[codon_table[0] == amino].iterrows():
        temp_triplets.append((row[1]).strip().upper())
        temp_values.append(row[2])
    
    return temp_triplets[temp_values.index(max(temp_values))]


def codon_optimize(sequence_id, cds, codon_table, force):
    """
    Checks if a given sequence starts with the start codon ('ATG') and promts you if it doesn't.
    Optimizes codons to most frequently used from the codon table.
    """
    run = True

    if cds[:3] != "ATG":
        if not force:
            run = prompt_force(sequence_id)
            force = True
    else:
        force = False

    if run:
        optimized_cds = ""

        for triplet in make_triplets(cds):
            
            if len(triplet) == 3:
                optimized_cds += reverse_translate(codon_table.loc[
                    codon_table[1] == triplet][0].iloc[0][0], codon_table)
            else:
                optimized_cds += triplet

        return optimized_cds, force
    
    else:
        return False, force


def translate_dna(sequence_id, cds, codon_table, force):
    """
    Translates coding DNA sequence to protein sequence.
    """
    run = True

    if cds[:3] != "ATG":
        if not force:
            run = prompt_force(sequence_id)
            force = True
    else:
        force = False

    if run:
        translation = ""

        for triplet in make_triplets(cds):
            
            if len(triplet) == 3:
                translation += codon_table.loc[codon_table[1] == triplet][0].iloc[0][0]
            else:
                translation += "?"

        return translation, force

    else:
        return False, force
