#cSpell: Disable#
import sys
import logging
import pyperclip
from random import choice
from seqtools.util import timeout
from seqtools.modules import Nucleotide, Protein
from seqtools.seqtools_config import LOGGING_CONFIG


@timeout(20)
def generate_dna(length, homopolymer=10, gc_stretch=20, restriction=False, ratio_gc=True,
                 protein=False):
    '''DNA generator; Generates random DNA sequence and executes basic tests
    TODO: gc_ratio has to be adjustable as is everything else, protein'''

    logger = logging.getLogger(__name__)

    dna = ('A', 'C', 'G', 'T')

    prot = ('A', 'C', 'D', 'E',
            'F', 'G', 'H', 'I',
            'K', 'L', 'M', 'N',
            'P', 'Q', 'R', 'S',
            'T', 'V', 'W', 'Y')

    restriction_enzymes = [
        "GGTCTC",  # BsaI
        "GAGACC",  # BsaI reverse
        "CGTCTC",  # BsmBI
        "GAGACG",  # BsmBI reverse
        "GCGGCCGC" # NotI
    ]

    SETTINGS = f'''
        length:            {length}
        homopolymer:    -n {homopolymer}
        gc_strech:      -g {gc_stretch}
        restriction:    -e {restriction}
        gc_ratio:       -r {ratio_gc}'''

    def generator(length, chars=dna):
        return ''.join(choice(chars) for _ in range(length))

    def check_restriction(sequence, restriction_set=restriction_enzymes):
        for restriction_site in restriction_set:
            if restriction_site in sequence:
                return True
        return False

    def check_homopolymer(sequence, upper_bound, chars=dna):
        for char in chars:
            if char * upper_bound in sequence:
                return True
        return False

    def check_gc_cont(sequence, upper_bound, lower_bound):
        if not upper_bound > sum(map(sequence.count, ('G', 'C'))) / len(sequence) > lower_bound:
            return True
        return False

    def check_gc_stretch(sequence, upper_bound):
        longest, at, gc = 0, 0, 0
        for char in sequence:
            if char in 'GC':
                gc += 1
            elif longest < gc:
                longest = gc
                gc = 0
            if char in 'AT':
                at += 1
            elif longest < at:
                longest = at
                at = 0
        if longest >= upper_bound:
            return True
        return False

    def run_tests(sequence_string, homopolymer=homopolymer, gc_stretch=gc_stretch,
                  restriction=restriction, ratio_gc=ratio_gc):
        a = check_homopolymer(sequence_string, upper_bound=homopolymer)
        b = check_gc_stretch(sequence_string, upper_bound=gc_stretch)
        c = False
        d = False
        if restriction:
            c = check_restriction(sequence_string)
        if ratio_gc:
            d = check_gc_cont(sequence_string, upper_bound=.6, lower_bound=.4)
        if a or b or c or d:
            return True
        return False

    try:
        if protein:
            seq = Protein('generated', f'M{generator(length-1, chars=prot)}*')
        else:
            candidates = []
            confirmed = []
            if length < 100:
                candidates.append(generator(length))
            else:
                while length > 100:
                    candidates.append(generator(100))
                    length -= 100
                candidates.append(generator(length))

            while candidates:
                if run_tests(candidates[0]):
                    candidates.append(generator(len(candidates[0])))
                    candidates.pop(0)
                else:
                    confirmed.append(candidates[0])
                    candidates.pop(0)    

            seq = Nucleotide('Generated DNA', ''.join(confirmed))

        logger.info(f'Generated sequence: {seq.sequence}')
        sys.stdout.write(seq.sequence)
        pyperclip.copy(seq.sequence)
        logger.info('Sequence copied to clipboard')
        return seq

    except TimeoutError:
        logger.error(f'''TIMEOUT ERROR! SETTINGS: {SETTINGS}''')
        print("Please try different settings.")

    except KeyboardInterrupt:
        logger.error(f'''CANCELED BY USER!''')
