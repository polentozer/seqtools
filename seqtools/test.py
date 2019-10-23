#cSpell: Disable#
import os
import re
import sys
import signal
import pandas
import random
import logging
import pyperclip
import logging.config
import matplotlib.pyplot as plt
from distutils.util import strtobool

### SEABORN ###
plt.style.use('seaborn-deep')
# plt.style.use('ggplot')

###### DATA ######
GGA_PART_TYPES = {
    'type1': {
        'prefix': 'GCATCGTCTCATCGGAGTCGGTCTCACCCT',
        'suffix': 'AACGAGAGACCAGCAGACCAGAGACGGCAT',
        'info': 'Left side assembly connector'
    },
    'type2': {
        'prefix': 'GCATCGTCTCATCGGTCTCAAACG',
        'suffix': 'TATGAGAGACCTGAGACGGCAT',
        'info': 'Promotor'
    },
    'type3': {
        'prefix': 'GCATCGTCTCATCGGTCTCAT',
        'suffix': 'ATCCAGAGACCTGAGACGGCAT',
        'info': 'CDS'
    },
    'type3a': {
        'prefix': 'GCATCGTCTCATCGGTCTCAT',
        'suffix': 'GGTTCTAGAGACCTGAGACGGCAT',
        'info': 'N-terminal CDS'
    },
    'type3b': {
        'prefix': 'GCATCGTCTCATCGGTCTCATTCT',
        'suffix': 'GGATCCAGAGACCTGAGACGGCAT',
        'info': 'CDS'
    },
    'type3t': {
        'prefix': 'GCATCGTCTCATCGGTCTCAT',
        'suffix': 'GGATCCTGAGACCTGAGACGGCAT',
        'info': 'True type3 CDS (GS linker, no STOP)'
    },
    'type4': {
        'prefix': 'GCATCGTCTCATCGGTCTCAATCC',
        'suffix': 'GCTGAGAGACCTGAGACGGCAT',
        'info': 'Terminator'
    },
    'type4a': {
        'prefix': 'GCATCGTCTCATCGGTCTCAATCC',
        'suffix': 'TGGCAGAGACCTGAGACGGCAT',
        'info': 'C-terminal CDS'
    },
    'type4b': {
        'prefix': 'GCATCGTCTCATCGGTCTCATGGC',
        'suffix': 'GCTGAGAGACCTGAGACGGCAT',
        'info': 'Terminator'
    },
    'type5': {
        'prefix': 'GCATCGTCTCATCGGAGTCGGTCTCAGCTG',
        'suffix': 'TACAAGAGACCAGCAGACCAGAGACGGCAT',
        'info': 'Right side assembly connector'
    },
    'type6': {
        'prefix': 'GCATCGTCTCATCGGTCTCATACA',
        'suffix': 'GAGTAGAGACCTGAGACGGCAT',
        'info': 'Yeast marker'
    },
    'type7': {
        'prefix': 'GCATCGTCTCATCGGTCTCAGAGT',
        'suffix': 'CCGAAGAGACCTGAGACGGCAT',
        'info': "3'-homology or yeast origin"
    },
    'type8': {
        'prefix': 'GCATCGTCTCATCGGTCTCACCGA',
        'suffix': 'CCCTAGAGACCAGAGACGGCAT',
        'info': 'E. coli marker and origin'
    },
    'type8a': {
        'prefix': 'GCATCGTCTCATCGGTCTCACCGA',
        'suffix': 'CAATAGAGACCAGAGACGGCAT',
        'info': 'E. coli marker and origin'
    },
    'type8b': {
        'prefix': 'GCATCGTCTCATCGGTCTCACAAT',
        'suffix': 'CCCTAGAGACCAGAGACGGCAT',
        'info': "5'-homology"
    },
    'typeX': {
        'prefix': 'GCATCGTCTCATCGGTCTCANNNN',
        'suffix': 'NNNNAGAGACCAGAGACGGCAT',
        'info': 'Custom parts'
    }
}

RESTRICTION_ENZYMES = {
    'BsaI': {
        'substrate': 'DNA',
        'recognition': 'GGTCTC',
        'jump': 1,
        'overhang': 4,
        'incubation_temperature': 37,
        'overhang_type': '5`',
        'methylation_sensitivity': {
            'Dam': False,
            'Dcm': True,
            'EcoKI': False
        },
        'description': 'Type 2 restriction enzyme used in modular cloning or MoClo for short. Sticky ends from different BsaI sites may not be compatible. BsaI can be used between 37 and 50 °C.'
    },
    'BsmBI': {
        'substrate': 'DNA',
        'recognition': 'CGTCTC',
        'jump': 1,
        'overhang': 4,
        'incubation_temperature': 55,
        'overhang_type': '5`',
        'methylation_sensitivity': {
            'Dam': False,
            'Dcm': False,
            'EcoKI': False
        },
        'description': 'Type 2 restriction enzyme used in modular cloning or MoClo for short. Sticky ends from different BsmBI sites may not be compatible.'
    },
    'NotI': {
        'substrate': 'DNA',
        'recognition': 'GCGGCCGC',
        'jump': 0,
        'overhang': 4,
        'incubation_temperature': 37,
        'overhang_type': '5`',
        'methylation_sensitivity': {
            'Dam': False,
            'Dcm': False,
            'EcoKI': False
        },
        'description': 'Classic restriction enzyme commonly used in cloning. This particular enzyme has has very rare recognition site.'
    },
    'BpiI': {
        'substrate': 'DNA',
        'recognition': 'GAAGAC',
        'jump': 2,
        'overhang': 4,
        'incubation_temperature': 37,
        'overhang_type': '5`',
        'methylation_sensitivity': {
            'Dam': False,
            'Dcm': False,
            'EcoKI': False
        },
        'description': 'Type 2 restriction enzyme used in modular cloning or MoClo for short. Sticky ends from different BsmBI sites may not be compatible.'
    },
    'EcoRI': {
        'substrate': 'DNA',
        'recognition': 'GAATTC',
        'jump': 0,
        'overhang': 4,
        'incubation_temperature': 37,
        'overhang_type': '5`',
        'methylation_sensitivity': {
            'Dam': False,
            'Dcm': False,
            'EcoKI': False
        },
        'description': 'Classic restriction enzyme commonly used in cloning.'
    }
}

LOGGING_CONFIG = {
    'version': 1,
    'formatters': {
        'simple': {
            'format': '%(asctime)s %(filename)s %(name)s %(levelname)s %(message)s'
        },
    },
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'level': 'INFO',
            'formatter': 'simple',
            'stream': 'ext://sys.stdout',
        },
        'file': {
            'class': 'logging.FileHandler',
            'level': 'DEBUG',
            'formatter': 'simple',
            'filename': 'seqtools.log'
        }
    },
    'loggers': {
        '__main__': {
            'level': 'DEBUG',
            'handlers': ['console', 'file'],
            'propagate': False
        }
    },
    'root': {
        'level': 'DEBUG',
        'handlers': ['console']
    }
}

COMMON_SPECIES = {
    'ecoli': '37762',
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
    'ATA', 'ATC', 'ATT', 'ATG', 'TGG', 'TAA', 'TAG', 'TGA'
]

STANDARD_GENETIC_CODE = [
    'R', 'R', 'R', 'R', 'R', 'R', 'L', 'L',
    'L', 'L', 'L', 'L', 'S', 'S', 'S', 'S',
    'S', 'S', 'T', 'T', 'T', 'T', 'P', 'P',
    'P', 'P', 'A', 'A', 'A', 'A', 'G', 'G',
    'G', 'G', 'V', 'V', 'V', 'V', 'K', 'K',
    'N', 'N', 'Q', 'Q', 'H', 'H', 'E', 'E',
    'D', 'D', 'Y', 'Y', 'C', 'C', 'F', 'F',
    'I', 'I', 'I', 'M', 'W', '*', '*', '*'
]

### LOGGER ###
logging.config.dictConfig(LOGGING_CONFIG)
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

###########################################
###                                     ###
###########################################
####                UTIL                ###
###########################################
###                                     ###
###########################################

####### TIMEOUT DECORATOR #######
class TimeoutError(Exception):
    def __init__(self, value='Timed Out'):
        self.value = value
    def __str__(self):
        return repr(self.value)

def timeout(seconds_before_timeout):
    def decorate(f):
        def handler(signum, frame):
            raise TimeoutError()
        def new_f(*args, **kwargs):
            old = signal.signal(signal.SIGALRM, handler)
            signal.alarm(seconds_before_timeout)
            try:
                result = f(*args, **kwargs)
            finally:
                # reinstall the old signal handler
                signal.signal(signal.SIGALRM, old)
                # cancel the alarm
                signal.alarm(0)
            return result
        new_f.__name__ = f.__name__
        return new_f
    return decorate
#################################


def bool_user_prompt(nucleotide_obj, operation):
    '''Prompts user for YES/NO response. Version 2.'''
    print(f'Sequence with ID {nucleotide_obj.sequence_id} is not a CDS. {operation} anyway? [y/n]')
    while True:
        try:
            if strtobool(input().lower()):
                nucleotide_obj.sequence_id += '|FORCED'
            return nucleotide_obj
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


def load_codon_table(species=None, taxonomy_id=None, custom=False, return_name=False):
    '''Load a codon table based on the organism's species ID'''
    # CODON_USAGE_DB = f'{os.path.dirname(__file__)}/data/codon_usage.spsum'
    # CUSTOM_CODON_USAGE_DB = f'{os.path.dirname(__file__)}/data/custom_table.spsum'
    CODON_USAGE_DB = './data/codon_usage.spsum'
    CUSTOM_CODON_USAGE_DB = './data/custom_table.spsum'

    if species in COMMON_SPECIES:
        taxonomy_id = COMMON_SPECIES[species]

    if custom:
        codon_usage = CUSTOM_CODON_USAGE_DB
    else:
        codon_usage = CODON_USAGE_DB

    with open(codon_usage) as cu:
        for header in cu:
            codon_counts = cu.readline()

            taxid, species = header.strip().split(':')[:2]

            if taxonomy_id:
                taxonomy_id = str(taxonomy_id)

            if taxonomy_id and taxonomy_id != taxid:
                continue

            table = list(zip(CODONS, STANDARD_GENETIC_CODE, [int(x) for x in codon_counts.split()]))
            table = pandas.DataFrame(table, columns=['Triplet', 'AA', 'Number'])
            table.set_index(['AA', 'Triplet'], inplace=True)
            table.sort_index(inplace=True)
            total = sum(table['Number'])

            table['Fraction'] = table.groupby('AA').transform(lambda x: x / x.sum())
            table['Frequency'] = table['Number'] / total * 1000

            break
    
    if return_name:
        return table, species

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
    x = codons.Fraction.cumsum() / codons.Fraction.cumsum().max() < random.random()

    return codons.iloc[x.sum()].name


def codon_table_10plus(table):
    '''Return a codon table only representing codons with > 10% occurrence frequency.'''

    table = table[table.Fraction >= 0.1]
    table = table.groupby(level=0).transform(lambda x: x / x.sum())

    return table

###########################################
###                FASTA                ###
###########################################

def open_fasta(file_input_paths, protein=False):
    '''Function for opening `.fasta` files. Files can contain multiple sequences.
    Returns list of sequence objects.'''
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


def write_fasta(sequence_list, path=''):
    '''Simple function for writing `.fasta` files.'''
    with open(f'{path}output.fasta', 'w') as file_out:
        for sequence in sequence_list:
            if sequence:
                file_out.write(f'\n{sequence}\n')


# NOTE: This will go away once code is refractured
DEFAULT_TABLE = load_codon_table(species='yali')

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
        return f'Sequence: >{self.sequence_id} {self.sequence}'

    def __str__(self):
        return self.sequence

    def __len__(self):
        return len(self.sequence)

    def __add__(self, other):
        return Sequence('concat', self.sequence + other.sequence)

    @property
    def fasta(self):
        return f'>{self.sequence_id}\n{self.sequence}'

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
        return Protein('concat', self.sequence + other.sequence)

    def __repr__(self):
        return f'Protein_Sequence: >{self.sequence_id} {self.sequence}'

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

    def reverse_translate(self, table=DEFAULT_TABLE, maximum=False):
        '''Returns optimized DNA sequence'''
        dna_sequence = list()

        if maximum:
            name = '|NUC-MAX'
            maximum = True
        else:
            name = '|NUC'
            maximum = False

        for amino in self.sequence:
            if amino == '?':
                dna_sequence.append('NNN')
            else:
                codons = table.loc[amino]
                dna_sequence.append(get_codon(codons, maximum=maximum))

        return Nucleotide(f'{self.sequence_id}{name}', ''.join(dna_sequence))


class Nucleotide(Sequence):
    '''NUCLEOTIDE sequence object'''

    def __init__(self, sequence_id, sequence):
        super().__init__(sequence_id, sequence)

    def __add__(self, other):
        return Nucleotide('concat', self.sequence + other.sequence)

    def __repr__(self):
        return f'Nucleotide_Sequence: >{self.sequence_id} {self.sequence}'

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
            return len(self) % 3 == 0

        def start(self):
            return self.sequence[:3] == 'ATG'

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

    def translate(self, table=DEFAULT_TABLE, check=False):
        '''Translate DNA sequence in PROTEIN sequence'''
        if not check:
            if not self.basic_cds:
                self = bool_user_prompt(self, 'Translate')
                if 'FORCED' not in self.sequence_id:
                    return self

        seq_id = self.sequence_id
        translation = list()
        table = table.reset_index(level='Triplet')

        for triplet in self.make_triplets():
            if len(triplet) == 3 and 'N' not in triplet:
                translation.append(table[table['Triplet'] == triplet].index[0])
            else:
                translation.append('?')

        return Protein(f'{seq_id}|PROT', ''.join(translation))

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

    def recode_sequence(self, replace, table=DEFAULT_TABLE, maximum=False):
        '''Recode a sequence to replace certain sequences using a given codon table.'''
        position = self.sequence.find(replace)

        if position < 0:
            return self

        position -= position % 3

        for i in range(position, position + (len(replace) // 3 + 1) * 3, 3):
            codon = self.sequence[i:i+3]
            options = table.loc[table.xs(codon, level=1).index[0]]

            if options.shape[0] == 1:
                continue

            if options.shape[0] > 0:
                new_codon = get_codon(options, maximum=maximum, recode=True, skip=[codon])
                break

        logger.warning(f'{codon} --> {new_codon}')

        if '|REC' not in self.sequence_id:
            self.sequence_id += '|REC'
        self.sequence = f'{self.sequence[:i]}{new_codon}{self.sequence[i+3:]}'

        return self

    def remove_cutsites(self, restriction_enzymes, table=DEFAULT_TABLE):
        '''Remove recognition sites for restriction enzymes.'''
        changes = 0

        for enzyme in restriction_enzymes:
            for cutsite in enzyme.cutsite_list:
                while cutsite in self.sequence:
                    logger.warning(f'{enzyme.enzyme_id} cuts ({cutsite})')
                    changes += 1
                    self = self.recode_sequence(cutsite, table=table)

        logger.info(f'remove_cutsites made {changes} changes in {self.sequence_id} sequence')

        return self

    def optimize_codon_usage(self, table=DEFAULT_TABLE, maximum=False):
        '''Optimize codon usage of a given DNA sequence'''
        if not self.basic_cds:
            self = bool_user_prompt(self, 'Optimize')
            if 'FORCED' not in self.sequence_id:
                return self

        seq_id = self.sequence_id
        optimized = self.translate(table=table).reverse_translate(table=table, maximum=maximum)

        return Nucleotide(f'{seq_id}|OPT', optimized.sequence)

    def make_part(self, part_type='3t', part_options=GGA_PART_TYPES):
        '''Make DNA part out of a given sequence'''
        logger.debug('Making parts...')
        seq_id = f'part_gge{part_type}_{self.sequence_id}'
        part = part_options[f'type{part_type}']
        if part_type in ('3t', '3a', '3b') and self.translate(check=True).sequence[-1] == '*':
            sequence = f'{part["prefix"]}{self.sequence[:-3]}{part["suffix"]}'
        else:
            sequence = f'{part["prefix"]}{self.sequence}{part["suffix"]}'

        return Nucleotide(seq_id, sequence)
    
    def harmonize(self, table_source, mode=0, table=DEFAULT_TABLE):
        '''Optimize codon usage of a given DNA sequence
        mode: 0 for closest frequency; 1 for same index'''
        if not self.basic_cds:
            self = bool_user_prompt(self, 'Special optimize')
            if 'FORCED' not in self.sequence_id:
                return self
        
        seq_id = self.sequence_id
        optimized = list()

        for amino, triplet in zip(self.translate(table=table).sequence, self.make_triplets()):
            if amino == '?':
                optimized.append('NNN')
            else:
                codons = table.loc[amino]
                source_codons = table_source.loc[amino]

                sorted_codons = sorted(codons['Fraction'])
                source_codon_freq = source_codons.loc[triplet]['Fraction']

                if mode == 0:
                    best, freq = 1, 0
                    for cod in sorted_codons:
                        current_best = abs(cod - source_codon_freq)
                        if current_best < best:
                            best, freq = current_best, cod

                    closest_freq_codon = codons[codons['Fraction'] == freq].index[0]

                    optimized.append(closest_freq_codon)
                
                elif mode == 1:
                    sorted_source_codons = sorted(source_codons['Fraction'])
                    source_codon_index = sorted_source_codons.index(source_codons.loc[triplet]['Fraction'])
                    same_index_codon = codons[codons['Fraction'] == sorted_codons[source_codon_index]].index[0]

                    optimized.append(same_index_codon)
                
                else:
                    print('ERROR')
                    return self
        
        return Nucleotide(f'{seq_id}|HAR{mode}', ''.join(optimized))
    
    def graph_codon_usage(self, window=16, other=None, other_id=None, table=DEFAULT_TABLE, minmax=True, target='Yarrowia lipolytica'):
        '''Graph codon frequency of a given gene'''

        if not self.basic_cds:
            return

        if isinstance(other, Nucleotide) and other.basic_cds:
            if other_id:
                if other_id in COMMON_SPECIES:
                    other_id = COMMON_SPECIES[other_id]
                table_other, species = load_codon_table(taxonomy_id=other_id, return_name=True)
            else:
                table_other = table
            if minmax:
                data = [x for x in zip(self.data_minmax(table=table, window=window), other.data_minmax(table=table_other, window=window))]
            else:
                data = [x for x in zip(self.data_fraction(table=table, window=window), other.data_fraction(table=table_other, window=window))]
        else:
            if minmax:
                data = self.data_minmax(table=table, window=window)
            else:
                data = self.data_fraction(table=table, window=window)

        x = range(len(data))
        zeros = [0 for i in x]

        if other:
            y1 = [i[0] for i in data]
            y2 = [i[1] for i in data]
            _, (ax0, ax1) = plt.subplots(2, 1, sharex=True, figsize=(12, 5))
            plt.subplots_adjust(left=.08, right=0.98, hspace=.5)

            ax0.plot(x, y1, alpha=0.8, linewidth=.5)
            ax0.set_title(f'Codon usage plot for {self.sequence_id} in {target}')

            if minmax:
                ax0.set_ylim(-100, 100)
                ax0.axhline(0, color='black', linewidth=.5)
                ax0.fill_between(x, y1, zeros, where=[True if y > 0 else False for y in y1], alpha=0.5, interpolate=True, color='C0')
                ax0.fill_between(x, y1, zeros, where=[True if y < 0 else False for y in y1], alpha=0.5, interpolate=True, color='C2')
                ax0.set_ylabel('%MinMax Value')
            else:
                ax0.set_ylabel('Fraction')

            ax1.plot(x, y2, alpha=0.8, linewidth=.5)
            if other_id:
                target = species
            ax1.set_title(f'Codon usage plot for {other.sequence_id} in {target}')

            if minmax:
                ax1.set_ylim(-100, 100)
                ax1.axhline(0, color='black', linewidth=.5)
                ax1.fill_between(x, y2, zeros, where=[True if y > 0 else False for y in y2], alpha=0.5, interpolate=True, color='C0')
                ax1.fill_between(x, y2, zeros, where=[True if y < 0 else False for y in y2], alpha=0.5, interpolate=True, color='C2')
                ax1.set_ylabel('%MinMax Value')
            else:
                ax1.set_ylabel('Fraction')

        else:
            _, ax = plt.subplots(1, 1, figsize=(12, 2))
            plt.subplots_adjust(left=.08, right=0.98, bottom=.25)
            ax.plot(x, data, alpha=0.8, linewidth=.5)
            ax.set_title(f'Codon usage plot for {self.sequence_id} in {target}')

            if minmax:
                ax.set_ylim(-100, 100)
                ax.axhline(0, color='black', linewidth=.5)
                ax.fill_between(x, data, zeros, where=[True if y > 0 else False for y in data], alpha=0.5, interpolate=True, color='C0')
                ax.fill_between(x, data, zeros, where=[True if y < 0 else False for y in data], alpha=0.5, interpolate=True, color='C2')
                ax.set_ylabel('%MinMax Value')
            else:
                ax.set_ylabel('Fraction')

        plt.xlim(-4, len(data)+4)
        plt.xlabel('Codon')
        plt.show()

        return

    def data_minmax(self, table=DEFAULT_TABLE, window=16):
        values, data = [], []
        tri_table = table.reset_index(level='Triplet')

        for triplet in self.make_triplets():
            freq = tri_table[tri_table['Triplet'] == triplet]['Frequency'][0]
            codons = table.loc[tri_table[tri_table['Triplet'] == triplet].index[0]]

            values.append((freq, max(codons.Frequency), min(codons.Frequency), sum(codons.Frequency)/len(codons)))

        for n in range(len(values)+1-window):
            temp = values[n:n+window]
            actual = sum([f[0] for f in temp]) / window
            maximum = sum([f[1] for f in temp]) / window
            minimum = sum([f[2] for f in temp]) / window
            average = sum([f[3] for f in temp]) / window

            maxi = ((actual - average) / (maximum - average)) * 100
            mini = ((average - actual) / (average - minimum)) * 100

            if maxi > 0:
                data.append(maxi)
            elif mini > 0:
                data.append(-mini)

        return data

    def data_fraction(self, table=DEFAULT_TABLE, window=16):
        values, data = [], []
        codons = table.reset_index().set_index(['Triplet'])

        for triplet in self.make_triplets():
            values.append(codons.loc[triplet]['Fraction'])

        for n in range(len(values)+1-window):
            data.append(sum([f for f in values[n:n+window]]) / window)

        return data


class Enzyme:
    def __init__(self, enzyme_id, substrate, enzyme_description=None):
        self.enzyme_id = enzyme_id
        self.substrate = substrate
        self.enzyme_description = enzyme_description

    def __repr__(self):
        return f'Enzyme: {self.enzyme_id}'

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
        return f'Restriction_Enzyme: {self.enzyme_id} ({self.recognition_sequence})'

    def __str__(self):
        return f'{self.enzyme_id} >> {self.recognition_sequence.sequence}' 

    @property
    def cutsite_list(self):
        return [self.recognition_sequence.sequence, self.recognition_sequence.reverse_complement.sequence]


###########################################
###                                     ###
###########################################
###             GENERATOR               ###
###########################################
###                                     ###
###########################################

@timeout(20)
def generate_dna(length, homopolymer=10, gc_stretch=20, restriction=False, ratio_gc=True,
                 protein=False):
    '''DNA generator; Generates random DNA sequence and executes basic tests
    TODO: gc_ratio has to be adjustable as is everything else, protein'''

    dna = ('A', 'C', 'G', 'T')

    prot = ('A', 'C', 'D', 'E',
            'F', 'G', 'H', 'I',
            'K', 'L', 'M', 'N',
            'P', 'Q', 'R', 'S',
            'T', 'V', 'W', 'Y')

    re = [
        "GGTCTC",  # BsaI
        "GAGACC",  # BsaI reverse
        "CGTCTC",  # BsmBI
        "GAGACG",  # BsmBI reverse
        "GCGGCCGC" # NotI
    ]

    setting = f'''
        length:            {length}
        homopolymer:    -n {homopolymer}
        gc_strech:      -g {gc_stretch}
        restriction:    -t {restriction}
        gc_ratio:       -r {ratio_gc}'''

    def generator(length, chars=dna):
        return ''.join(random.choice(chars) for _ in range(length))

    def check_restriction(sequence, restriction_set=re):
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
        logger.error(f'''TIMEOUT ERROR! SETTINGS: {setting}''')
        print("Please try different settings.")

    except KeyboardInterrupt:
        logger.error(f'''CANCELED BY USER!''')

###########################################
###                                     ###
###########################################
###            TEST SUITE               ###
###########################################
###                                     ###
###########################################

if __name__ == '__main__':
    '''testing stuff'''

    loc_new_table = './data/codon_usage.spsum'

    test_dna = Nucleotide('dna', 'ATGGCCCGAAAGGCCCCACACATCGACTAA')
    test_dna2 = Nucleotide('idi1', 'ATGACTGCCGACAACAATAGTATGCCCCATGGTGCAGTATCTAGTTACGCCAAATTAGTGCAAAACCAAACACCTGAAGACATTTTGGAAGAGTTTCCTGAAATTATTCCATTACAACAAAGACCTAATACCCGATCTAGTGAGACGTCAAATGACGAAAGCGGAGAAACATGTTTTTCTGGTCATGATGAGGAGCAAATTAAGTTAATGAATGAAAATTGTATTGTTTTGGATTGGGACGATAATGCTATTGGTGCCGGTACCAAGAAAGTTTGTCATTTAATGGAAAATATTGAAAAGGGTTTACTACATCGTGCATTCTCCGTCTTTATTTTCAATGAACAAGGTGAATTACTTTTACAACAAAGAGCCACTGAAAAAATAACTTTCCCTGATCTTTGGACTAACACATGCTGCTCTCATCCACTATGTATTGATGACGAATTAGGTTTGAAGGGTAAGCTAGACGATAAGATTAAGGGCGCTATTACTGCGGCGGTGAGAAAACTAGATCATGAATTAGGTATTCCAGAAGATGAAACTAAGACAAGGGGTAAGTTTCACTTTTTAAACAGAATCCATTACATGGCACCAAGCAATGAACCATGGGGTGAACATGAAATTGATTACATCCTATTTTATAAGATCAACGCTAAAGAAAACTTGACTGTCAACCCAAACGTCAATGAAGTTAGAGACTTCAAATGGGTTTCACCAAATGATTTGAAAACTATGTTTGCTGACCCAAGTTACAAGTTTACGCCTTGGTTTAAGATTATTTGCGAGAATTACTTATTCAACTGGTGGGAGCAATTAGATGACCTTTCTGAAGTGGAAAATGACAGGCAAATTCATAGAATGCTATAA')
    test_prot = Protein('prot', 'MARKAPHID*')

    ecmad = Nucleotide('EcMAD', 'ATGAACAACGGCACAAATAATTTTCAGAACTTCATCGGGATCTCAAGTTTGCAGAAAACGCTGCGCAATGCTCTGATCCCCACGGAAACCACGCAACAGTTCATCGTCAAGAACGGAATAATTAAAGAAGATGAGTTACGTGGCGAGAACCGCCAGATTCTGAAAGATATCATGGATGACTACTACCGCGGATTCATCTCTGAGACTCTGAGTTCTATTGATGACATAGATTGGACTAGCCTGTTCGAAAAAATGGAAATTCAGCTGAAAAATGGTGATAATAAAGATACCTTAATTAAGGAACAGACAGAGTATCGGAAAGCAATCCATAAAAAATTTGCGAACGACGATCGGTTTAAGAACATGTTTAGCGCCAAACTGATTAGTGACATATTACCTGAATTTGTCATCCACAACAATAATTATTCGGCATCAGAGAAAGAGGAAAAAACCCAGGTGATAAAATTGTTTTCGCGCTTTGCGACTAGCTTTAAAGATTACTTCAAGAACCGTGCAAATTGCTTTTCAGCGGACGATATTTCATCAAGCAGCTGCCATCGCATCGTCAACGACAATGCAGAGATATTCTTTTCAAATGCGCTGGTCTACCGCCGGATCGTAAAATCGCTGAGCAATGACGATATCAACAAAATTTCGGGCGATATGAAAGATTCATTAAAAGAAATGAGTCTGGAAGAAATATATTCTTACGAGAAGTATGGGGAATTTATTACCCAGGAAGGCATTAGCTTCTATAATGATATCTGTGGGAAAGTGAATTCTTTTATGAACCTGTATTGTCAGAAAAATAAAGAAAACAAAAATTTATACAAACTTCAGAAACTTCACAAACAGATTCTATGCATTGCGGACACTAGCTATGAGGTCCCGTATAAATTTGAAAGTGACGAGGAAGTGTACCAATCAGTTAACGGCTTCCTTGATAACATTAGCAGCAAACATATAGTCGAAAGATTACGCAAAATCGGCGATAACTATAACGGCTACAACCTGGATAAAATTTATATCGTGTCCAAATTTTACGAGAGCGTTAGCCAAAAAACCTACCGCGACTGGGAAACAATTAATACCGCCCTCGAAATTCATTACAATAATATCTTGCCGGGTAACGGTAAAAGTAAAGCCGACAAAGTAAAAAAAGCGGTTAAGAATGATTTACAGAAATCCATCACCGAAATAAATGAACTAGTGTCAAACTATAAGCTGTGCAGTGACGACAACATCAAAGCGGAGACTTATATACATGAGATTAGCCATATCTTGAATAACTTTGAAGCACAGGAATTGAAATACAATCCGGAAATTCACCTAGTTGAATCCGAGCTCAAAGCGAGTGAGCTTAAAAACGTGCTGGACGTGATCATGAATGCGTTTCATTGGTGTTCGGTTTTTATGACTGAGGAACTTGTTGATAAAGACAACAATTTTTATGCGGAACTGGAGGAGATTTACGATGAAATTTATCCAGTAATTAGTCTGTACAACCTGGTTCGTAACTACGTTACCCAGAAACCGTACAGCACGAAAAAGATTAAATTGAACTTTGGAATACCGACGTTAGCAGACGGTTGGTCAAAGTCCAAAGAGTATTCTAATAACGCTATCATACTGATGCGCGACAATCTGTATTATCTGGGCATCTTTAATGCGAAGAATAAACCGGACAAGAAGATTATCGAGGGTAATACGTCAGAAAATAAGGGTGACTACAAAAAGATGATTTATAATTTGCTCCCGGGTCCCAACAAAATGATCCCGAAAGTTTTCTTGAGCAGCAAGACGGGGGTGGAAACGTATAAACCGAGCGCCTATATCCTAGAGGGGTATAAACAGAATAAACATATCAAGTCTTCAAAAGACTTTGATATCACTTTCTGTCATGATCTGATCGACTACTTCAAAAACTGTATTGCAATTCATCCCGAGTGGAAAAACTTCGGTTTTGATTTTAGCGACACCAGTACTTATGAAGACATTTCCGGGTTTTATCGTGAGGTAGAGTTACAAGGTTACAAGATTGATTGGACATACATTAGCGAAAAAGACATTGATCTGCTGCAGGAAAAAGGTCAACTGTATCTGTTCCAGATATATAACAAAGATTTTTCGAAAAAATCAACCGGGAATGACAACCTTCACACCATGTACCTGAAAAATCTTTTCTCAGAAGAAAATCTTAAGGATATCGTCCTGAAACTTAACGGCGAAGCGGAAATCTTCTTCAGGAAGAGCAGCATAAAGAACCCAATCATTCATAAAAAAGGCTCGATTTTAGTCAACCGTACCTACGAAGCAGAAGAAAAAGACCAGTTTGGCAACATTCAAATTGTGCGTAAAAATATTCCGGAAAACATTTATCAGGAGCTGTACAAATACTTCAACGATAAAAGCGACAAAGAGCTGTCTGATGAAGCAGCCAAACTGAAGAATGTAGTGGGACACCACGAGGCAGCGACGAATATAGTCAAGGACTATCGCTACACGTATGATAAATACTTCCTTCATATGCCTATTACGATCAATTTCAAAGCCAATAAAACGGGTTTTATTAATGATAGGATCTTACAGTATATCGCTAAAGAAAAAGACTTACATGTGATCGGCATTGATCGGGGCGAGCGTAACCTGATCTACGTGTCCGTGATTGATACTTGTGGTAATATAGTTGAACAGAAAAGCTTTAACATTGTAAACGGCTACGACTATCAGATAAAACTGAAACAACAGGAGGGCGCTAGACAGATTGCGCGGAAAGAATGGAAAGAAATTGGTAAAATTAAAGAGATCAAAGAGGGCTACCTGAGCTTAGTAATCCACGAGATCTCTAAAATGGTAATCAAATACAATGCAATTATAGCGATGGAGGATTTGTCTTATGGTTTTAAAAAAGGGCGCTTTAAGGTCGAACGGCAAGTTTACCAGAAATTTGAAACCATGCTCATCAATAAACTCAACTATCTGGTATTTAAAGATATTTCGATTACCGAGAATGGCGGTCTCCTGAAAGGTTATCAGCTGACATACATTCCTGATAAACTTAAAAACGTGGGTCATCAGTGCGGCTGCATTTTTTATGTGCCTGCTGCATACACGAGCAAAATTGATCCGACCACCGGCTTTGTGAATATCTTTAAATTTAAAGACCTGACAGTGGACGCAAAACGTGAATTCATTAAAAAATTTGACTCAATTCGTTATGACAGTGAAAAAAATCTGTTCTGCTTTACATTTGACTACAATAACTTTATTACGCAAAACACGGTCATGAGCAAATCATCGTGGAGTGTGTATACATACGGCGTGCGCATCAAACGTCGCTTTGTGAACGGCCGCTTCTCAAACGAAAGTGATACCATTGACATAACCAAAGATATGGAGAAAACGTTGGAAATGACGGACATTAACTGGCGCGATGGCCACGATCTTCGTCAAGACATTATAGATTATGAAATTGTTCAGCACATATTCGAAATTTTCCGTTTAACAGTGCAAATGCGTAACTCCTTGTCTGAACTGGAGGACCGTGATTACGATCGTCTCATTTCACCTGTACTGAACGAAAATAACATTTTTTATGACAGCGCGAAAGCGGGGGATGCACTTCCTAAGGATGCCGATGCAAATGGTGCGTATTGTATTGCATTAAAAGGGTTATATGAAATTAAACAAATTACCGAAAATTGGAAAGAAGATGGTAAATTTTCGCGCGATAAACTCAAAATCAGCAATAAAGATTGGTTCGACTTTATCCAGAATAAGCGCTATCTCTAA')

    new_table = load_codon_table(species='yali')

    re_list = []
    for renzyme in RESTRICTION_ENZYMES.keys():
        re_list.append(
            Restriction_Enzyme(renzyme,
                RESTRICTION_ENZYMES[renzyme]['substrate'],
                RESTRICTION_ENZYMES[renzyme]['recognition'],
                jump=RESTRICTION_ENZYMES[renzyme]['jump'],
                overhang=RESTRICTION_ENZYMES[renzyme]['overhang'],
                enzyme_description=RESTRICTION_ENZYMES[renzyme]['description']
            )
        )

    print(re_list)



    # print(test_dna2)
    # for _ in range(10):
    #     print(test_dna2.optimize_codon_usage(maximum=False).remove_cutsites(re_list))
    # print(test_dna2.remove_cutsites(re_list))

    # print(test_dna2.optimize_codon_usage(maximum=False).remove_cutsites(re_list).make_part('3t').fasta)

    print()
    print()

    # ecmad.optimize_codon_usage().graph_codon_usage(other=ecmad)
    # l1 = [1, 2, 3, 4]
    # l2 = [2, 4, 6, 8]

    # x = list(map(lambda x: x*2, (l1, l2)))

    # print(x)

    # Salvia officinalis: 38868 table
    # s_table = load_codon_table(taxonomy_id=38868)

    # print(test_dna2.harmonize(s_table, mode=1))

    # print(test_dna2.translate())
    # for x in test_dna2.translate().sequence:
    #     print(x)

    # for amino, cod in zip(test_dna2.translate().sequence, test_dna2.make_triplets()):
    #     print(cod, amino)

    # print(new_table)
    # print(s_table)


    ########################################################################
    # amino = 'V'
    # src = 'GTA'
    # codons = new_table.loc[amino]
    # scodons = s_table.loc[amino]

    # print(codons)
    # print(scodons)

    # print()

    # for x, y in codons['Fraction'].items():
    #     print(x, y)

    # print()

    # sorted_codons = sorted(codons['Fraction'])
    # sorted_scodons = sorted(scodons['Fraction'])

    # scodon_occ = scodons.loc[src]['Fraction']
    # scodon_index = sorted_scodons.index(scodon_occ)
    
    # maped = map(lambda x, y: [x, y], sorted_codons, sorted_scodons)

    # # for x in maped:
    # #     print(x)
    # print(list(maped))
    # print(scodon_occ)
    # print(scodon_index)

    # print(sorted_codons)
    # print(sorted_scodons)
    # # print(sorted_codons[-1])
    # # print(codons[codons['Fraction'] == sorted_codons[-1]])

    # print("\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n")

    # best, occ = 1, 0
    # for cod in sorted_codons:
    #     current = abs(cod - scodon_occ)
    #     print(f'best: {best:.5f}    occ: {occ:.5f}    current: {current:.5f}    cod: {cod:.5f}')
    #     if current < best:
    #         best, occ = current, cod
    
    # closest_freq = codons[codons['Fraction'] == occ]
    # same_index = codons[codons['Fraction'] == sorted_codons[scodon_index]]

    # print(best)
    # print(occ)

    # print(f'Closest frequency codon: {closest_freq.index[0]}')
    # print(f'Same index codon: {same_index.index[0]}')
    ########################################################################

    source_table = load_codon_table(species='ecoli')
    # yeast_table = load_codon_table(species='yeast')

    # options = new_table.reset_index().set_index(['Triplet'])
    # source_options = source_table.reset_index().set_index(['Triplet'])

    # test_dna2_opt = test_dna2.optimize_codon_usage(maximum=True)
    # test_dna2_opt = test_dna2.optimize_codon_usage(maximum=False)
    # test_dna2_opt = test_dna2.harmonize(table_source=source_table)

    # mad = Nucleotide('mad', 'ATGAACAACGGCACAAATAATTTTCAGAACTTCATCGGGATCTCAAGTTTGCAGAAAACGCTGCGCAATGCTCTGATCCCCACGGAAACCACGCAACAGTTCATCGTCAAGAACGGAATAATTAAAGAAGATGAGTTACGTGGCGAGAACCGCCAGATTCTGAAAGATATCATGGATGACTACTACCGCGGATTCATCTCTGAGACTCTGAGTTCTATTGATGACATAGATTGGACTAGCCTGTTCGAAAAAATGGAAATTCAGCTGAAAAATGGTGATAATAAAGATACCTTAATTAAGGAACAGACAGAGTATCGGAAAGCAATCCATAAAAAATTTGCGAACGACGATCGGTTTAAGAACATGTTTAGCGCCAAACTGATTAGTGACATATTACCTGAATTTGTCATCCACAACAATAATTATTCGGCATCAGAGAAAGAGGAAAAAACCCAGGTGATAAAATTGTTTTCGCGCTTTGCGACTAGCTTTAAAGATTACTTCAAGAACCGTGCAAATTGCTTTTCAGCGGACGATATTTCATCAAGCAGCTGCCATCGCATCGTCAACGACAATGCAGAGATATTCTTTTCAAATGCGCTGGTCTACCGCCGGATCGTAAAATCGCTGAGCAATGACGATATCAACAAAATTTCGGGCGATATGAAAGATTCATTAAAAGAAATGAGTCTGGAAGAAATATATTCTTACGAGAAGTATGGGGAATTTATTACCCAGGAAGGCATTAGCTTCTATAATGATATCTGTGGGAAAGTCAATTCTTTTATGAACCTGTATTGTCAGAAAAATAAAGAAAACAAAAATTTATACAAACTTCAGAAACTTCACAAACAGATTCTATGCATTGCGGACACTAGCTATGAGGTCCCGTATAAATTTGAAAGTGACGAGGAAGTGTACCAATCAGTTAACGGCTTCCTTGATAACATTAGCAGCAAACATATAGTCGAAAGATTACGCAAAATCGGCGATAACTATAACGGCTACAACCTGGATAAAATTTATATCGTGTCCAAATTTTACGAGAGCGTTAGCCAAAAAACCTACCGCGACTGGGAAACAATTAATACCGCCCTCGAAATTCATTACAATAATATCTTGCCGGGTAACGGTAAAAGTAAAGCCGACAAAGTAAAAAAAGCGGTTAAGAATGATTTACAGAAATCCATCACCGAAATAAATGAACTAGTGTCAAACTATAAGCTGTGCAGTGACGACAACATCAAAGCGGAGACTTATATACATGAGATTAGCCATATCTTGAATAACTTTGAAGCACAGGAATTGAAATACAATCCGGAAATTCACCTAGTTGAATCCGAGCTCAAAGCGAGTGAGCTTAAAAACGTGCTGGACGTGATCATGAATGCGTTTCATTGGTGTTCGGTTTTTATGACTGAGGAACTTGTTGATAAAGACAACAATTTTTATGCGGAACTGGAGGAGATTTACGATGAAATTTATCCAGTAATTAGTCTGTACAACCTGGTTCGTAACTACGTTACCCAGAAACCGTACAGCACGAAAAAGATTAAATTGAACTTTGGAATACCGACGTTAGCAGACGGTTGGTCAAAGTCCAAAGAGTATTCTAATAACGCTATCATACTGATGCGCGACAATCTGTATTATCTGGGCATCTTTAATGCGAAGAATAAACCGGACAAGAAGATTATCGAGGGTAATACGTCAGAAAATAAGGGTGACTACAAAAAGATGATTTATAATTTGCTCCCGGGTCCCAACAAAATGATCCCGAAAGTTTTCTTGAGCAGCAAGACGGGGGTGGAAACGTATAAACCGAGCGCCTATATCCTAGAGGGGTATAAACAGAATAAACATATCAAATCTTCAAAAGACTTTGATATCACTTTCTGTCATGATCTGATCGACTACTTCAAAAACTGTATTGCAATTCATCCCGAGTGGAAAAACTTCGGTTTTGATTTTAGCGACACCAGTACTTATGAGGACATTTCCGGGTTTTATCGTGAGGTAGAGTTACAAGGTTACAAGATTGATTGGACATACATTAGCGAAAAAGACATTGATCTGCTGCAGGAAAAAGGTCAACTGTATCTGTTCCAGATATATAACAAAGATTTTTCGAAAAAATCAACCGGGAATGACAACCTTCACACCATGTACCTGAAAAATCTTTTCTCAGAAGAAAATCTTAAGGATATCGTCCTGAAACTTAACGGCGAAGCGGAAATCTTCTTCAGGAAGAGCAGCATAAAGAACCCAATCATTCATAAAAAAGGCTCGATTTTAGTCAACCGTACCTACGAAGCAGAAGAAAAAGACCAGTTTGGCAACATTCAAATTGTGCGTAAAAATATTCCGGAAAACATTTATCAGGAGCTGTACAAATACTTCAACGATAAAAGCGACAAAGAGCTGTCTGATGAAGCAGCCAAACTGAAGAATGTAGTGGGACACCACGAGGCAGCGACGAATATAGTCAAGGACTATCGCTACACGTATGATAAATACTTCCTTCATATGCCTATTACGATCAATTTCAAAGCCAATAAAACGGGTTTTATTAATGATAGGATCTTACAGTATATCGCTAAAGAAAAAGACTTACATGTGATCGGCATTGATCGGGGCGAGCGTAACCTGATCTACGTGTCCGTGATTGATACTTGTGGTAATATAGTTGAACAGAAAAGCTTTAACATTGTAAACGGCTACGACTATCAGATAAAACTGAAACAACAGGAGGGCGCTAGACAGATTGCGCGGAAAGAATGGAAAGAAATTGGTAAAATTAAAGAGATCAAAGAGGGCTACCTGAGCTTAGTAATCCACGAGATCTCTAAAATGGTAATCAAATACAATGCAATTATAGCGATGGAGGATTTGTCTTATGGTTTTAAAAAAGGGCGCTTTAAGGTCGAACGGCAAGTTTACCAGAAATTTGAAACCATGCTCATCAATAAACTCAACTATCTGGTATTTAAAGATATTTCGATTACCGAGAATGGCGGACTCCTGAAAGGTTATCAGCTGACATACATTCCTGATAAACTTAAAAACGTGGGTCATCAGTGCGGCTGCATTTTTTATGTGCCTGCTGCATACACGAGCAAAATTGATCCGACCACCGGCTTTGTGAATATCTTTAAATTTAAAGACCTGACAGTGGACGCAAAACGTGAGTTCATTAAAAAATTTGACTCAATTCGTTATGACAGTGAAAAAAATCTGTTCTGCTTTACATTTGACTACAATAACTTTATTACGCAAAACACGGTCATGAGCAAATCATCGTGGAGTGTGTATACATACGGCGTGCGCATCAAACGTCGCTTTGTGAACGGCCGCTTCTCAAACGAAAGTGATACCATTGACATAACCAAAGATATGGAGAAAACGTTGGAAATGACGGACATTAACTGGCGCGATGGCCACGATCTTCGTCAAGACATTATAGATTATGAAATTGTTCAGCACATATTCGAAATTTTCCGTTTAACAGTGCAAATGCGTAACTCCTTGTCTGAACTGGAGGACCGTGATTACGATAGACTCATTTCACCTGTACTGAACGAAAATAACATTTTTTATGACAGCGCGAAAGCGGGGGATGCACTTCCTAAGGATGCCGATGCAAATGGTGCGTATTGTATTGCATTAAAAGGGTTATATGAAATTAAACAAATTACCGAAAATTGGAAAGAAGATGGTAAATTTTCGCGCGATAAACTCAAAATCAGCAATAAAGATTGGTTCGACTTTATCCAGAATAAGCGCTATCTCTCTCGAGCCGACCCCAAGAAGAAGCGAAAGGTGTAA')
    # ylmad = ecmad.optimize_codon_usage(table=DEFAULT_TABLE)
    # ylmadm = ecmad.optimize_codon_usage(maximum=True)
    # scmad = ecmad.optimize_codon_usage(table=yeast_table)
    # ylmad1 = ecmad.harmonize(table_source=source_table)
    # scmad1 = ecmad.harmonize(table=yeast_table, table_source=source_table)
    # comp = []
    # mad.graph_codon_usage()
    # scmad1.graph_codon_usage()
    # ylmad.graph_codon_usage(scmad)
    # scmad.graph_codon_usage(ylmad)
    # ylmadm.graph_codon_usage(ecmad)
    # ecmad.graph_codon_usage(ylmad, table=source_table)
    # ecmad.graph_codon_usage(ylmad)
    # ylmad1.graph_harmonization(ecmad, other_id='ecoli')
    # ylmad1.graph_codon_usage(other=ecmad, other_id='ecoli')
    # ylmad1.graph_codon_usage(other=ecmad)
    # ylmad1.graph_codon_usage()
    # ylmad.graph_codon_usage()
    # ecmad.graph_codon_usage()

    # print(sum(yeast_table['Number']))
    # print(source_table)

    # ty = DEFAULT_TABLE.reset_index(level='Triplet')

    # data = []

    # for triplet in test_dna2.make_triplets():
    #     row = ty[ty['Triplet'] == triplet]
    #     amino = row.index[0]
    #     freq = row['Frequency'][0]

    #     codons = DEFAULT_TABLE.loc[amino]

    #     data.append((freq, max(codons.Frequency), min(codons.Frequency), sum(codons.Frequency)/len(codons)))

    # new_data = []
    # window = 16
    # for n in range(len(data)+1-window):
    #     temp = data[n:n+window]
    #     actual = sum([f[0] for f in temp]) / window
    #     maximum = sum([f[1] for f in temp]) / window
    #     minimum = sum([f[2] for f in temp]) / window
    #     average = sum([f[3] for f in temp]) / window

    #     maxi = ((actual - average) / (maximum - average)) * 100
    #     mini = ((average - actual) / (average - minimum)) * 100

    #     # print(f'maxi: {maxi:.1f} | mini: {mini:.1f}')

    #     if maxi > 0:
    #         new_data.append(maxi)
    #     elif mini > 0:
    #         new_data.append(-mini)

    # nnnn = pandas.DataFrame(test_dna2.minmax(), columns=['MinMax'])
    
    # print(nnnn)
    # nnnn.plot()
    # plt.show()



    # # for source_triplet, new_triplet in zip(test_dna2.make_triplets(), test_dna2_opt.make_triplets()):
    # for source_triplet, new_triplet in zip(ecmad.make_triplets(), ylmad1.make_triplets()):

    #     codon_freq = options.loc[new_triplet]['Fraction']
    #     source_codon_freq = options.loc[source_triplet]['Fraction']
    #     native_codon_freq = source_options.loc[source_triplet]['Fraction']
    #     comp.append((source_codon_freq, codon_freq, native_codon_freq))
    #     # print(new_triplet, codon_freq, '|', source_triplet, source_codon_freq)

    # c2 = []
    # temp = []
    # while comp:
    #     freq = comp[0]
    #     if (len(temp) == 6 or len(comp) < 6) and temp:
    #         source_window = sum([x[0] for x in temp])/len(temp)
    #         new_window = sum([x[1] for x in temp])/len(temp)
    #         native_window = sum([x[2] for x in temp])/len(temp)
    #         c2.append((source_window, new_window, native_window))
    #         temp = []
    #     temp.append(freq)
    #     comp.pop(0)

    # # tab = pandas.DataFrame(comp, columns=['source', 'new', 'native'])
    # tab = pandas.DataFrame(c2, columns=['source', 'new', 'native'])

    # print(tab)

    # opt_value = sum(tab['new']) - sum(tab['source'])

    # print(opt_value)

    # # tab.plot(y='source')
    # # tab.plot(y='new')
    # # tab.plot(y=['source', 'new'])
    # tab.plot(y=['native', 'new'])
    # plt.show()

    # r = 20
    # q = 4

    # t = [x for x in range(r)]

    # for a in range(r-q+1):
    #     b = t[a:a+q]
    #     print(b)

    # for k in GGA_PART_TYPES.keys():
    #     print(k[4:])
    
    if '3a' in [k[4:] for k in GGA_PART_TYPES.keys()]:
        print('True')
    
    x = '''>sfGFP
ATGAGAAAGGGTGAGGAGCTGTTTACAGGCGTTGTGCCCATCCTTGTGGAGCTTGACGG
CGATGTCAACGGTCATAAGTTCTCCGTTCGAGGAGAAGGCGAGGGTGATGCCACCAACG
GAAAGCTCACCCTTAAGTTCATCTGCACTACCGGTAAGCTGCCCGTGCCTTGGCCTACC
CTGGTCACCACTCTCACCTATGGTGTGCAGTGCTTCGCCCGATACCCCGATCACATGAA
GCAGCATGACTTTTTTAAGTCTGCCATGCCTGAAGGCTACGTCCAGGAGCGAACCATCT
CTTTCAAGGATGATGGAACCTACAAGACCCGAGCTGAGGTCAAGTTTGAGGGTGATACT
CTGGTCAACCGAATCGAGCTGAAGGGAATTGATTTTAAGGAGGATGGAAACATCCTCGG
TCACAAGCTGGAGTACAACTTCAACTCCCACAACGTCTACATCACCGCTGATAAGCAGA
AGAACGGAATTAAGGCTAACTTCAAGATTCGACACAATGTCGAAGACGGTTCCGTTCAG
CTCGCTGACCACTACCAGCAGAACACTCCTATTGGTGACGGACCCGTCCTGCTGCCCGA
CAACCATTACCTCTCGACCCAGTCGGTCCTTTCTAAGGACCCTAACGAGAAGCGAGATC
ACATGGTCCTTCTGGAGTTCGTCACTGCCGCTGGAATTACCCACGGCATGGATGAGCTT
TACAAG'''

for line in x.splitlines():
    print(line)
