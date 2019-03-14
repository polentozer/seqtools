def make_triplets(string):
    """
    Makes chunks of 3 characters from a long string.
    """
    return [string[start:start + 3] for start in range(0, len(string), 3)]


def dna_operation(sequences, codon_table, force, optimize=None, analyze=None):
    """
    Optimizes DNA triplets to use most common codons from the codon_table.
    """
    solutions = {}

    if analyze:
        for sequence_id, sequence in sequences.items():
            optimization_value(sequence_id, sequence, codon_table, force)

    elif optimize:
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
    Optimizes codons to most frequently used from the codon table.
    """
    run, force = check_start(cds[:3], sequence_id, force)

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
        return None, force


def translate_dna(sequence_id, cds, codon_table, force):
    """
    Translates coding DNA sequence to protein sequence.
    """
    run, force = check_start(cds[:3], sequence_id, force)

    if run:
        translation = ""

        for triplet in make_triplets(cds):

            if len(triplet) == 3:
                translation += codon_table.loc[codon_table[1] == triplet][0].iloc[0][0]
            else:
                translation += "?"

        return translation, force

    else:
        return None, force


def check_start(first_codon, sequence_id, force):
    """
    Checks if first codon equals to 'ATG' and prompts to force operation, if it doesn't.
    """
    if first_codon != "ATG":
        if not force:
            cont_prompt = input(
                "Sequence with ID `{0}` is not a CDS, optimize anyway? (Y/n): ".format(sequence_id))
            if cont_prompt.upper() == "Y" or cont_prompt.upper() == "YES" or cont_prompt == "":
                run, force = True, True
            else:
                run, force = False, False
        else:
            run, force = True, True
    else:
        run, force = True, False

    return run, force


def optimization_value(sequence_id, cds, codon_table, force):
    """
    Calculates the codon optimization value of a given sequence.
    MAX VALUE: 1
    MIN VALUE: 0
    """
    run, force = check_start(cds[:3], sequence_id, force)

    if run:
        values = []

        for amino, original, optimized in zip(
                translate_dna(sequence_id, cds, codon_table, force)[0], make_triplets(cds),
                make_triplets(codon_optimize(sequence_id, cds, codon_table, force)[0])):

            if original == optimized:
                x = 1
            else:
                x = 0
            values.append(x)
            print(amino, original, optimized, x)

        if len(values) == sum(values):
            print(f"{sequence_id} is fully optimized!")

        else:
            opt_val = sum(values) / len(values)
            print(f"Optimization percentage for {sequence_id} is {opt_val*100:.2f} %")

        print()

    else:
        return None, force


def remove_type2_sites(sequence_id, cds, codon_table):
    """
    Finds and removes type2 restriction sites

    BsaI:
    GGTCTCN^NNNN N
    CCAGAGN NNNN^N

    BsmBI:
    CGTCTCN^NNNN N
    GCAGAGN NNNN^N
    """
    pass
