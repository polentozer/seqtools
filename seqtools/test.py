#cSpell: Disable#
import os
import re
import sys
import pandas
import logging
import pyperclip
import numpy as np
from random import choice
from collections import deque
from distutils.util import strtobool


### LOGGER ###
logger = logging.getLogger(__name__)
# logger.setLevel(logging.DEBUG)
# logging.basicConfig(filename='seqtools.log', level=logging.DEBUG)
logging.basicConfig(level=logging.DEBUG)


###########################################
###                                     ###
###########################################
####                UTIL                ###
###########################################
###                                     ###
###########################################

def codon_table_parser(path_to_table):
    '''Used for loading old style codon usage table
    TODO: deprecated'''
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


def melting_temperature(dna_sequence):
    '''Calculate and return the Tm using the "Wallace rule".

    Tm = 4 degC * (G+C) + 2 degC * (A+T)

    The Wallace rule (Thein & Wallace 1986, in Human genetic diseases: a
    practical approach, 33-50) is often used as rule of thumb for approximate
    Tm calculations for primers of 14 to 20 nt length.

    Non-dNA characters (e.g. E, F, J, !, 1, etc) are ignored in this method.
    '''
    weak = ('A', 'T', 'W')
    strong = ('C', 'G', 'S')
    return 2 * sum(map(dna_sequence.count, weak)) + 4 * sum(map(dna_sequence.count, strong))


def load_codon_table(species=None, taxonomy_id=None, custom=False):
    '''Load a codon table based on the organism's species ID'''
    # CODON_USAGE_DB = f'{os.path.dirname(__file__)}/data/codon_usage.spsum'
    # CUSTOM_CODON_USAGE_DB = f'{os.path.dirname(__file__)}/data/custom_table.spsum'
    CODON_USAGE_DB = './data/codon_usage.spsum'
    CUSTOM_CODON_USAGE_DB = './data/custom_table.spsum'

    COMMON_SPECIES = {
        'ecoli': '83333',
        'yeast':  '4932',
        'human': '9606',
        'bsub': '1432',
        'yali': '284591'
    }

    CODONS = [
        'CGA', 'CGC', 'CGG', 'CGT', 'AGA', 'AGG', 'CTA', 'CTC',
        'CTG', 'CTT', 'TTA', 'TTG', 'TCA', 'TCC', 'TCG', 'TCT',
        'AGC', 'AGT', 'ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'CCC',
        'CCG', 'CCT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC',
        'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'AAA', 'AAG',
        'AAC', 'AAT', 'CAA', 'CAG', 'CAC', 'CAT', 'GAA', 'GAG',
        'GAC', 'GAT', 'TAC', 'TAT', 'TGC', 'TGT', 'TTC', 'TTT',
        'ATA', 'ATC', 'ATT', 'ATG', 'TGG', 'TAA', 'TAG', 'TGA']

    STANDARD_GENETIC_CODE = [
        'R', 'R', 'R', 'R', 'R', 'R', 'L', 'L',
        'L', 'L', 'L', 'L', 'S', 'S', 'S', 'S',
        'S', 'S', 'T', 'T', 'T', 'T', 'P', 'P',
        'P', 'P', 'A', 'A', 'A', 'A', 'G', 'G',
        'G', 'G', 'V', 'V', 'V', 'V', 'K', 'K',
        'N', 'N', 'Q', 'Q', 'H', 'H', 'E', 'E',
        'D', 'D', 'Y', 'Y', 'C', 'C', 'F', 'F',
        'I', 'I', 'I', 'M', 'W', '*', '*', '*']

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

            table = list(
                zip(CODONS, STANDARD_GENETIC_CODE, [int(x) for x in codon_counts.split()]))
            table = pandas.DataFrame(table, columns=['Triplet', 'AA', 'Number'])
            table.set_index(['AA', 'Triplet'], inplace=True)
            table.sort_index(inplace=True)

            table['Fraction'] = table.groupby('AA').transform(lambda x: x / x.sum())
            break
    
    return table


def get_codon(codons, maximum=False, recode=False, skip=[]):
    '''Returns a "locally-optimized" codon. Locally-optimized = mimics the
    codon frequency in the table. Maximum uses the most common codon.'''

    if recode:
        codons = codons.loc[[cod for cod in codons.index if cod not in skip]]

    # Deterministic allocation of codon based on the highest frequency
    if maximum:
        return codons.Fraction.idxmax()
    
    # Stochastic allocation of codon
    x = codons.Fraction.cumsum() / codons.Fraction.cumsum().max() < np.random.rand()

    return codons.iloc[x.sum()].name


def codon_table_10plus(table):
    """Return a codon table only representing codons with > 10% occurrence frequency."""

    table = table[table.Fraction >= 0.1]
    table = table.groupby(level=0).transform(lambda x: x / x.sum())

    return table

###########################################
###                FASTA                ###
###########################################

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

# NOTE: This will go away once code is refractured
default_table = load_codon_table(species='yali')

###########################################
###                                     ###
###########################################
###             MODULES                 ###
###########################################
###                                     ###
###########################################    


class Sequence:
    '''Biological sequence object'''

    def __init__(self, sequence_id, sequence):
        self.sequence_id = sequence_id
        self.sequence = sequence.upper()

    def __repr__(self):
        return f'>{self.sequence_id}\n{self.sequence}'

    def __str__(self):
        return self.sequence

    def __len__(self):
        return len(self.sequence)

    def kmer_occurrence(self, threshold, length=8):
        kmers = {}

        for i in range(len(self) - length + 1):
            kmer = self.sequence[i:i + length]
            if kmer not in kmers:
                kmers[kmer] = 0
            kmers[kmer] += 1

        return [a for a in sorted(kmers.items(), key=lambda x: x[1]) if a[1] > threshold][::-1]


class Protein(Sequence):
    '''PROTEIN sequence object'''

    def __init__(self, sequence_id, sequence):
        super().__init__(sequence_id, sequence)

    def __add__(self, other):
        return Protein('concatenated', self.sequence + other.sequence)

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, string):
        allowed_characters = re.compile(r'[^\*\?GALMFWKQESPVICYHRNDT]')
        if sequence_match(string, allowed_characters.search):
            self._sequence = string
        else:
            raise ValueError

    def reverse_translate(self, table=default_table, maximum=False):
        '''Returns optimized DNA sequence'''
        dna_sequence = list()

        if maximum:
            name = '|NUC-MAX'
            maximum = True
        else:
            name = '|NUC'
            maximum = False

        for amino in self.sequence:
            codons = table.loc[amino]
            dna_sequence.append(get_codon(codons, maximum=maximum))

        return Nucleotide(f'{self.sequence_id}{name}', ''.join(dna_sequence))


class Nucleotide(Sequence):
    '''NUCLEOTIDE sequence object'''

    def __init__(self, sequence_id, sequence):
        super().__init__(sequence_id, sequence)

    def __add__(self, other):
        return Nucleotide('concatenated', self.sequence + other.sequence)

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, string):
        allowed_characters = re.compile(r'[^ACTGNU]')
        if not sequence_match(string, allowed_characters.search):
            raise ValueError
        self._sequence = string

    @property
    def basic_cds(self):
        '''Returns True if sequence is CDS or false if its not'''
        if self.sequence[:3] == 'ATG' and len(self) % 3 == 0:
            return True
        else:
            return False

    def check_cds(self):
        '''Checks CDS'''
        def triplet(self):
            return len(self.sequence) % 3 == 0
        
        def start(self):
            return self.sequence[:3] == "ATG"
        
        def stop(self):
            prot = self.translate(check=True)
            return prot.sequence[-1] == "*"
        
        def no_internal_stop(self):
            prot = self.translate(check=True)
            return not "*" in prot.sequence[:-1]
        
        tests = [triplet, start, stop, no_internal_stop]
        result = True
        for test in tests:
            if not test(self):
                logger.warning(f'{self.sequence_id} failed on {test.__name__}')
                result = False
        
        return result

    @property
    def reverse_complement(self):
        '''Returns reverse complement of given DNA sequence'''
        return Nucleotide(f'{self.sequence_id}|REVC',
            self.sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1])

    def make_triplets(self):
        '''Makes list of chunks 3 characters long from a sequence'''
        return [self.sequence[start:start + 3] for start in range(0, len(self.sequence), 3)]
    
    def translate(self, table=default_table, check=False):
        '''Translate DNA sequence in PROTEIN sequence'''
        if not check:
            prompt = f'Sequence with ID {self.sequence_id} is not a CDS. Translate anyway?'

            if not self.basic_cds:
                if not bool_user_prompt(prompt):
                    return None
                else:
                    new_sequence_id = f'{self.sequence_id}|FORCED'
            else:
                new_sequence_id = self.sequence_id
        else:
            new_sequence_id = self.sequence_id

        translation = list()
        table = table.reset_index(level='Triplet')

        for triplet in self.make_triplets():
            if len(triplet) == 3:
                translation.append(table[table['Triplet'] == triplet].index[0])
            else:
                translation.append('?')

        return Protein(f'{new_sequence_id}|PROT', ''.join(translation))
    
    # @property
    # def translation(self):
    #     '''Returns protein sequence for easier workflow
    #     TODO: I don't know if i wanna keep this or not'''
    #     return self.translate().sequence

    def check_common_errors(self):
        logger.info(f'Checking common errors in sequence with ID "{self.sequence_id}"...')
        for nuc in 'ACGT':
            logger.info(f'Checking "{nuc}" homopolymer')
            if nuc*11 in self.sequence:
                logger.warning(f'ID "{self.sequence_id}" failed on long "{nuc}" homopolymer')
            else:
                logger.info("pass")
        logger.info(f'Checking for unknown base pairs')
        if 'N' in self.sequence:
            logger.warning(f'ID "{self.sequence_id}" failed on "N" (unknown base pair) in sequence')
        else:
                logger.info("pass")
        logger.info(f'Checking BsaI restriction sites')
        if 'GGTCTC' in self.sequence or 'GAGACC' in self.sequence:
            logger.warning(f'ID "{self.sequence_id}" failed on BsaI restriction sites')
        else:
                logger.info("pass")
        logger.info(f'Checking BsmBI restriction sites')
        if 'CGTCTC' in self.sequence or 'GAGACG' in self.sequence:
            logger.warning(f'ID "{self.sequence_id}" failed on BsmBI restriction sites')
        else:
                logger.info("pass")
        logger.info(f'Checking too short sequence')
        if len(self) < 10:
            logger.warning(f'ID "{self.sequence_id}" failed on small sequence (<10 bp)')
        else:
                logger.info("pass")
        logger.info(f'Checking too long sequence')
        if len(self) > 20000:
            logger.warning(f'ID "{self.sequence_id}" failed on long sequence (>20 kbp)')
        else:
                logger.info("pass")
        return

    def recode_sequence(self, replace, table=default_table, maximum=False):
        '''Recode a sequence to replace certain sequences using a given codon table.'''
        position = self.sequence.find(replace)

        if position < 0:
            return self
        
        position -= position % 3

        for i in range(position, position + (len(replace) // 3 + 1) * 3, 3):
            codon = self.sequence[i:i+3]
            options = table.loc[table.xs(codon, level=1).index[0]]

            if options.shape[0] > 0:
                new_codon = get_codon(options, maximum=maximum, recode=True, skip=[codon])
                break
        
        logger.warning(f'{codon} --> {new_codon}')

        if '|REC' in self.sequence_id:
            name = self.sequence_id
        else:
            name = f'{self.sequence_id}|REC'

        return Nucleotide(name, f'{self.sequence[:i]}{new_codon}{self.sequence[i+3:]}')

    # def remove_cutsites(self, cut_sites):
    #     '''Remove recognition sites for restriction enzymes.'''
    #     changes = 0

    #     for enzyme, seq in cut_sites + []

    def optimize_codon_usage(self, codon_table):
        '''Optimize codon usage given with codon usage table
        TODO: fix table (optimize_codon_usage)'''
        prompt = f'Sequence with ID {self.sequence_id} is not a CDS. Optimize anyway?'

        if not self.basic_cds:
            if not bool_user_prompt(prompt):
                return self
            else:
                new_sequence_id = f'{self.sequence_id}|FORCED'
        else:
            new_sequence_id = self.sequence_id

        optimized_sequence = ''
        for triplet in self.make_triplets():
            if len(triplet) == 3:
                temp_values = {}
                temp_amino = codon_table[codon_table[1] == triplet].iloc[0][0]

                for _, row in codon_table[codon_table[0] == temp_amino].iterrows():
                    temp_values[row[1]] = row[2]

                optimized_sequence += max(temp_values, key=temp_values.get)
            else:
                optimized_sequence += triplet

        return Nucleotide(f'{new_sequence_id}|OPT', optimized_sequence)

    def optimization_value(self, codon_table):
        '''
        Calculates the codon optimization value of a given sequence
        MAX VALUE: 1
        MIN VALUE: 0
        TODO: fix table (optimization_value)
        '''
        if not self.basic_cds:
            return 0

        values = []

        for _, original, optimized in zip(
                self.translate(codon_table).sequence,
                self.make_triplets(),
                self.optimize_codon_usage(codon_table).make_triplets()):

            x = 0
            if original == optimized:
                x = 1
            values.append(x)

        return 100 * sum(values) / len(values)


class Enzyme:
    def __init__(self, enzyme_id, substrate, enzyme_description=None):
        self.enzyme_id = enzyme_id
        self.substrate = substrate
        self.enzyme_description = enzyme_description
    
    def __repr__(self):
        return f'Enzyme: {self.enzyme_id}\n{self.enzyme_description}'
    
    def __str__(self):
        return self.enzyme_id


class Restriction_Enzyme(Enzyme):
    def __init__(self, enzyme_id, substrate, recognition_sequence, jump=0, overhang=0,
                 enzyme_description=None):
        super().__init__(enzyme_id, substrate, enzyme_description)
        self.recognition_sequence = Nucleotide(enzyme_id, recognition_sequence)
        self.jump = jump
        self.overhang = overhang

    def __repr__(self):
        if self.jump > 0:
            return f'Enzyme: {self.enzyme_id}\n{self.recognition_sequence}{"N"*self.jump}.{"N"*self.overhang}.N\n{self.enzyme_description}'
        else:
            return f'Enzyme: {self.enzyme_id}\n{self.recognition_sequence}\n\n\
                {self.enzyme_description}'

    def __str__(self):
        return f'{self.enzyme_id} >> {self.recognition_sequence.sequence}' 


###########################################
###                                     ###
###########################################
###             GENERATOR               ###
###########################################
###                                     ###
###########################################

def generate_dna(length, homopolymer=4, gc_stretch=6, type2=False, ratio_gc=False, protein=False):
    '''DNA generator
    Generates random DNA sequence and executes basic tests'''

    dna = ('A', 'C', 'G', 'T')
    prot = (
        'A', 'C', 'D', 'E',
        'F', 'G', 'H', 'I',
        'K', 'L', 'M', 'N',
        'P', 'Q', 'R', 'S',
        'T', 'V', 'W', 'Y')

    restriction_enzymes = [
        "GGTCTC",  # BsaI
        "GAGACC",  # BsaI reverse
        "CGTCTC",  # BsmBI
        "GAGACG",  # BsmBI reverse
        "TCTAGA",  # XbaI
        "GAATTC",  # EcoRI
    ]

    try:
        if not protein:
            while True:
                candidate = generator(length, chars=dna)

                a = check_homopolymer(candidate, upper_bound=homopolymer, chars=dna)
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
                    return dna
        else:
            return Protein('generated', f'M{generator(length-1, chars=prot)}*')
    
    except KeyboardInterrupt:
        print(f"""\nCanceled.\n
        SETTINGS:
            length:            {length}
            homopolymer:    -n {homopolymer}
            gc_strech:      -g {gc_stretch}
            w/o_type2s:     -t {type2}
            gc_ratio:       -r {ratio_gc}""")


def generator(length, chars):
    return ''.join(np.random.choice(chars) for _ in range(length))


def check_type2(sequence, restriction_set):
    for restriction_site in restriction_set:
        if restriction_site in sequence:
            return True
    return False


def check_homopolymer(sequence, upper_bound, chars):
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



if __name__ == '__main__':
    '''testing stuff'''

    restriction_enzymes = {
        'BsaI': {
            'substrate': 'DNA',
            'recognition': 'GGTCTC',
            'jump': 1,
            'overhang': 4,
            'description': 'Type 2 restriction enzyme used in modular cloning or MoClo for short.'
        },
        'BsmBI': {
            'substrate': 'DNA',
            'recognition': 'CGTCTC',
            'jump': 1,
            'overhang': 4,
            'description': 'Type 2 restriction enzyme used in modular cloning or MoClo for short.'
        }
    }

    gga_part_types = {
        'type1': {
            'prefix': 'GCATCGTCTCATCGGAGTCGGTCTCNCCCT',
            'suffix': 'AACGNGAGACCAGCAGACCAGAGACGGCAT',
            'info': 'Left side assembly connector'
        },
        'type2': {
            'prefix': 'GCATCGTCTCATCGGTCTCNAACG',
            'suffix': 'TATGNGAGACCTGAGACGGCAT',
            'info': 'Promotor'
        },
        'type3': {
            'prefix': 'GCATCGTCTCATCGGTCTCNT',
            'suffix': 'ATCCNGAGACCTGAGACGGCAT',
            'info': 'CDS'
        },
        'type3a': {
            'prefix': 'GCATCGTCTCATCGGTCTCNT',
            'suffix': 'GGTTCTNGAGACCTGAGACGGCAT',
            'info': 'N-terminal CDS'
        },
        'type3b': {
            'prefix': 'GCATCGTCTCATCGGTCTCNTTCT',
            'suffix': 'GGATCCNGAGACCTGAGACGGCAT',
            'info': 'CDS'
        },
        'type3t': {
            'prefix': 'GCATCGTCTCATCGGTCTCNT',
            'suffix': 'GGATCCTGAGACCTGAGACGGCAT',
            'info': 'True type3 CDS (GS linker, no STOP)'
        },
        'type4': {
            'prefix': 'GCATCGTCTCATCGGTCTCNATCC',
            'suffix': 'GCTGNGAGACCTGAGACGGCAT',
            'info': 'Terminator'
        },
        'type4a': {
            'prefix': 'GCATCGTCTCATCGGTCTCNATCC',
            'suffix': 'TGGCNGAGACCTGAGACGGCAT',
            'info': 'C-terminal CDS'
        },
        'type4b': {
            'prefix': 'GCATCGTCTCATCGGTCTCNTGGC',
            'suffix': 'GCTGNGAGACCTGAGACGGCAT',
            'info': 'Terminator'
        },
        'type5': {
            'prefix': 'GCATCGTCTCATCGGAGTCGGTCTCNGCTG',
            'suffix': 'TACANGAGACCAGCAGACCAGAGACGGCAT',
            'info': 'Right side assembly connector'
        },
        'type6': {
            'prefix': 'GCATCGTCTCATCGGTCTCNTACA',
            'suffix': 'GAGTNGAGACCTGAGACGGCAT',
            'info': 'Yeast marker'
        },
        'type7': {
            'prefix': 'GCATCGTCTCATCGGTCTCNGAGT',
            'suffix': 'CCGANGAGACCTGAGACGGCAT',
            'info': "3'-homology or yeast origin"
        },
        'type8': {
            'prefix': 'GCATCGTCTCATCGGTCTCNCCGA',
            'suffix': 'CCCTNGAGACCAGAGACGGCAT',
            'info': 'E. coli marker and origin'
        },
        'type8a': {
            'prefix': 'GCATCGTCTCATCGGTCTCNCCGA',
            'suffix': 'CAATNGAGACCAGAGACGGCAT',
            'info': 'E. coli marker and origin'
        },
        'type8b': {
            'prefix': 'GCATCGTCTCATCGGTCTCNCAAT',
            'suffix': 'CCCTNGAGACCAGAGACGGCAT',
            'info': "5'-homology"
        },
        'typeX': {
            'prefix': 'GCATCGTCTCATCGGTCTCNNNNN',
            'suffix': 'NNNNNGAGACCAGAGACGGCAT',
            'info': 'Custom parts'
        }
    }


    loc_old_table = './data/sample_table.csv'
    loc_new_table = './data/codon_usage.spsum'

    test_dna = Nucleotide('dna', 'ATGGCCCGAAAGGCCCCACACATCGACTAA')
    test_prot = Protein('prot', 'MARKAPHID*')
    
    old_table = codon_table_parser(loc_old_table)
    new_table = load_codon_table(species='yali')

    # print(test_dna)
    # print(test_prot)
    # print(repr(test_dna))
    # print(repr(test_prot))
    # print(len(test_dna))
    # print(len(test_prot))
    # print(test_dna.basic_cds)
    # print(test_dna.make_triplets())
    # print(test_dna.reverse_complement())
    # print(test_dna.reverse_complement().make_triplets())
    # print(test_dna.translate(old_table))
    # print(test_dna.reverse_complement().translate(old_table))
    # print(test_dna.reverse_complement().translate1(new_table))
    # print(test_dna.optimization_value(old_table))
    # print(test_dna.reverse_complement().optimization_value(old_table))
    # print(test_dna.kmer_occurrence(1, length=3))
    # print(test_dna.optimize_codon_usage(old_table))

    # print(new_table)
    # v_table = new_table.loc['V']
    # print(v_table)
    # print()
    # print(v_table.Fraction.cumsum())
    # print(v_table.Fraction.cumsum().max())
    # x = np.random.rand()
    # print(x)
    # print(
    #     v_table.iloc[
    #         (v_table.Fraction.cumsum() / v_table.Fraction.cumsum().max() < x).sum()
    #     ].name
    # )

    # freq = {
    #     "GTA": 0,
    #     "GTC": 0,
    #     "GTG": 0,
    #     "GTT": 0
    # }

    # for _ in range(10000):
    #     x = np.random.rand()
    #     chosen = v_table.iloc[(v_table.Fraction.cumsum() / v_table.Fraction.cumsum().max() < x).sum()].name
    #     freq[chosen] += 1
    
    # print(freq)
    # print()
    # max_triplet = v_table[v_table.Fraction == v_table.Fraction.max()].index[0]
    # max_triplet = v_table[v_table.index != 'GTG'].Fraction.idxmax()
    # max_triplet = v_table.Fraction.idxmax()

    # print(max_triplet)
    # print()

    # for _ in range(20):
        # print(test_prot.reverse_translate(new_table, maximum=True).sequence)
        # print(test_prot.reverse_translate(new_table, maximum=False).sequence)

    # x = new_table.index.Triplet('TAA')


    # x = new_table.reset_index(level='Triplet')

    # print(x)
    # print(x[x['Triplet'] == 'TAA'].index[0])

    # print(test_dna.translate(old_table))
    # print(test_dna.translate(new_table))

    # tt = codon_table_10plus(new_table)

    # print(tt)

    # print()
    # print()
    # codon = test_dna.sequence[3:6]
    # codon = 'GTG'
    # skip = ['GTA', 'GTT']
    # print(codon)

    # choices = codon_table_10plus(new_table).ix[new_table.ix[codon]]
    # choices = new_table.loc[new_table.xs(codon, level=1).index[0]]

    # table = new_table.reset_index(level='Triplet')

    # newcodon = 
    # print(options)
    # print(choices)
    # print()
    # print(list(choices.index))
    # print(type(choices.index))
    # print(choices[choices.index != skip[0]])
    # print(choices.loc[[x for x in choices.index if x not in skip]])
    # print()
    # print()
    # print()
    # print(test_dna)
    # print(test_dna.recode_sequence('GCCCCA', new_table))

    # print(test_dna.translate())

    # print(test_dna.check_cds())
    # test_dna.check_common_errors()
    # print()
    # print()
    # print()
    l = list(test_dna.sequence)
    d = deque(test_dna.sequence[1:])

    print(d)
    print(type(d))

    match = []
    longest = []

    while d:
        for i, item in enumerate(d):
            if l[i] == item:
                match.append(item)
            else:
                if len(longest) < len(match):
                    longest = match
                match = []
        d.popleft()
    
    print(longest)
    print(melting_temperature(test_dna.sequence))
    print()
    print(test_dna)
    print(repr(test_dna))
    print(test_prot)
    print(repr(test_prot))
    print()
    print()

    bsai = Restriction_Enzyme(
        'BsaI',
        restriction_enzymes['BsaI']['substrate'],
        restriction_enzymes['BsaI']['recognition'],
        jump=restriction_enzymes['BsaI']['jump'],
        overhang=restriction_enzymes['BsaI']['overhang'],
        enzyme_description=restriction_enzymes['BsaI']['description']
    )

    none_prot = Protein('none', '')

    print(repr(none_prot))

    print(bsai)
    print(bsai.jump)
    print(bsai.enzyme_description)
    print()
    print(repr(bsai))
