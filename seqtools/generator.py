import sys
import pyperclip
from random import choice
from seqtools.modules import Nucleotide


def generate_dna(length, single_repeats=4, gc_stretch=6, type2=False, ratio_gc=False):
    '''
    DNA generator
    '''
    dna = ['A', 'C', 'G', 'T']
    restriction_enzymes = [
        "GGTCTC",  # BsaI
        "GAGACC",  # BsaI reverse
        "CGTCTC",  # BsmBI
        "GAGACG",  # BsmBI reverse
        "TCTAGA",  # XbaI
        "GAATTC",  # EcoRI
    ]

    try:
        while True:
            candidate = generator(length, chars=dna)

            a = check_char_repeats(candidate, upper_bound=single_repeats, chars=dna)
            b = check_gc_stretch(candidate, upper_bound=gc_stretch)
            c = False
            d = False

            if type2:
                c = check_type2(candidate, restriction_set=restriction_enzymes)

            if ratio_gc:
                d = gc_cont(candidate, upper_bound=.6, lower_bound=.4)

            if not (a or b or c or d):
                dna = Nucleotide('Generated DNA', candidate)
                sys.stdout.write(dna.sequence)
                print("\n\nSequence has been copied to clipboard!")
                pyperclip.copy(dna.sequence)
                return
    
    except KeyboardInterrupt:
        print(f"""\nCanceled.\n
        SETTINGS:
            length:         {length}
            single_repeats: -n {single_repeats}
            gc_strech:      -g {gc_stretch}
            w/o_type2s:     -t {type2}
            gc_ratio:       -r {ratio_gc}""")


def generator(length, chars):
    result = ''
    while len(result) < length:
        result += choice(chars)
    return result


def check_type2(sequence, restriction_set):
    for restriction_site in restriction_set:
        if restriction_site in sequence:
            return True
    return False


def check_char_repeats(sequence, upper_bound, chars):
    for char in chars:
        if char * upper_bound in sequence:
            return True
    return False


def gc_cont(sequence, upper_bound, lower_bound):
    gc = 0
    for char in sequence:
        if char == "G" or char == "C":
            gc += 1

    ratio = gc / len(sequence)

    if not upper_bound > ratio > lower_bound:
        return True
    return False


def check_gc_stretch(sequence, upper_bound):
    longest = 0
    gc = 0
    at = 0
    for char in sequence:
        if char == "G" or char == "C":
            gc += 1
        else:
            if longest < gc:
                longest = gc
            gc = 0

        if char == "A" or char == "T":
            at += 1
        else:
            if longest < at:
                longest = at
            at = 0

    if longest >= upper_bound:
        return True
    return False
