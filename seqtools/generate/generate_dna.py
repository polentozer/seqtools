import pyperclip
from seqtools.modules import Nucleotide
from random import choice


def generate_dna(length, r=4, s=6, t=False, g=False):
    run = True
    while run:
        fin = generator(length, chars=dna)

        a = check_char_repeats(s=fin, upper_bound=r, chars=dna)
        d = check_gc_streach(s=fin, upper_bound=s)

        if not t:
            c = check_type2(s=fin, restriction_set=restriction_enzymes)
        else:
            c = False

        if not g:
            b = gc_cont(s=fin, upper_bound=.6, lower_bound=.4)
        else:
            g = False

        if a or b or c or d:
            run = True
        else:
            run = False

    else:
        print(fin)
        print("Sequence has been copied to clipboard!")
        pyperclip.copy(fin)


###
# UTIL
#

dna = ['A', 'C', 'G', 'T']
restriction_enzymes = [
    "GGTCTC",  # BsaI
    "GAGACC",  # BsaI reverse
    "CGTCTC",  # BsmBI
    "GAGACG",  # BsmBI reverse
    "TCTAGA",  # XbaI
    "GAATTC",  # EcoRI
]


def generator(l, chars):
    result = ''
    for i in range(l):
        result += choice(chars)

    return result


def check_char_repeats(s, upper_bound, chars):
    for c in chars:
        if c * upper_bound in s:
            return True
        else:
            return False


def gc_cont(s, upper_bound, lower_bound):
    gc = 0
    for c in s:
        if c == "G" or c == "C":
            gc += 1

    r = gc / len(s)

    if upper_bound > r > lower_bound:
        return False
    else:
        return True


def check_gc_streach(s, upper_bound):
    longest = 0
    gc = 0
    at = 0
    for c in s:
        if c == "G" or c == "C":
            gc += 1
        else:
            if longest < gc:
                longest = gc
            gc = 0

        if c == "A" or c == "T":
            at += 1
        else:
            if longest < at:
                longest = at
            at = 0

    if longest >= upper_bound:
        return True
    else:
        return False


def check_type2(s, restriction_set):
    for r in restriction_set:
        if r in s:
            return True
        else:
            return False
