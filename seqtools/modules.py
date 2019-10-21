#cSpell: Disable#
import re
import pandas
import logging
import matplotlib.pyplot as plt
from seqtools.seqtools_config import GGA_PART_TYPES, LOGGING_CONFIG, COMMON_SPECIES
from seqtools.util import bool_user_prompt, sequence_match, get_codon, load_codon_table

## LEGACY ##
DEFAULT_TABLE = load_codon_table(species='yali')
plt.style.use('seaborn-deep')


class Sequence:
    '''Biological sequence object'''

    def __init__(self, sequence_id, sequence, logger=None):
        self.sequence_id = sequence_id
        self.sequence = sequence.upper()
        self.logger = logger or logging.getLogger(__name__)

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
        return f'>{self.sequence_id}\n{self.sequence}\n'

    def kmer_analysis(self, threshold, length=8):
        self.logger.debug('Calculating k-mer occurrence...')
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
        allowed_characters = re.compile(r'[^\*\?GALMFWKQESPVICYHRNDTX]')
        if not sequence_match(string, allowed_characters.search):
            self.logger.error('Protein sequence includes unallowed character(s)')
            raise ValueError
        self._sequence = string

    def reverse_translate(self, table=DEFAULT_TABLE, maximum=False):
        '''Returns optimized DNA sequence'''
        self.logger.debug('Making reverse translation...')
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
        allowed_characters = re.compile(r'[^ACTGNUSW]')
        if not sequence_match(string, allowed_characters.search):
            self.logger.error('Nucleotide sequence includes unallowed character(s)')
            raise ValueError
        self._sequence = string

    @property
    def basic_cds(self):
        '''Returns True if sequence is CDS or false if its not'''
        self.logger.debug('Checking basic CDS...')
        if self.sequence[:3] == 'ATG' and len(self) % 3 == 0:
            return True
        self.logger.error(f'{self.sequence_id} is not a CDS')
        return False

    def check_cds(self):
        '''Checks CDS'''
        self.logger.debug('Checking CDS...')
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
                self.logger.error(f'{self.sequence_id} failed on {test.__name__}')
                result = False

        return result

    @property
    def reverse_complement(self):
        '''Returns reverse complement of given DNA sequence'''
        self.logger.debug('Making reverse complement...')
        return Nucleotide(f'{self.sequence_id}|REVC',
            self.sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1])

    def make_triplets(self):
        '''Makes list of chunks 3 characters long from a sequence'''
        self.logger.debug('Making triplets...')
        return [self.sequence[start:start+3] for start in range(0, len(self.sequence), 3)]

    def melting_temperature(self):
        '''Calculate and return the Tm using the "Wallace rule".

        Tm = 4°C * (G+C) + 2°C * (A+T)

        The Wallace rule (Thein & Wallace 1986, in Human genetic diseases: a
        practical approach, 33-50) is often used as rule of thumb for approximate
        Tm calculations for primers of 14 to 20 nt length.

        Non-dNA characters (e.g. E, F, J, !, 1, etc) are ignored in this method.
        '''
        self.logger.debug('Calculating melting temperature...')
        weak = ('A', 'T', 'W')
        strong = ('C', 'G', 'S')
        return 2*sum(map(self.sequence.count, weak)) + 4*sum(map(self.sequence.count, strong))

    def translate(self, table=DEFAULT_TABLE, check=False):
        '''Translate DNA sequence in PROTEIN sequence'''
        self.logger.debug('Making translation...')
        if not check:
            if not self.basic_cds:
                self = bool_user_prompt(self, 'Translate')
                if 'FORCED' not in self.sequence_id:
                    return self
                else:
                    self.logger.debug('Continuing with "FORCED"...')
        seq_id = self.sequence_id
        translation = list()
        table = table.reset_index(level='Triplet')
        for triplet in self.make_triplets():
            if len(triplet) == 3 and 'N' not in triplet:
                translation.append(table[table['Triplet'] == triplet].index[0])
            else:
                self.logger.warning(f'Unknown translation for codon: {triplet}')
                translation.append('?')

        return Protein(f'{seq_id}|PROT', ''.join(translation))

    def check_common_errors(self):
        self.logger.debug('Checking common errors on DNA sequence')
        self.logger.info(f'Checking common errors in sequence with ID "{self.sequence_id}"...')
        for nucleotide in 'ACGT':
            self.logger.info(f'Checking "{nucleotide}" homopolymer')
            if nucleotide * 11 in self.sequence:
                self.logger.warning(f'ID "{self.sequence_id}" failed on long "{nucleotide}" homopolymer')
            else:
                self.logger.info("pass")
        self.logger.info(f'Checking for unknown base pairs')
        if 'N' in self.sequence:
            self.logger.warning(f'ID "{self.sequence_id}" failed on "N" (unknown base pair) in sequence')
        else:
                self.logger.info("pass")
        self.logger.info(f'Checking BsaI restriction sites')
        if 'GGTCTC' in self.sequence or 'GAGACC' in self.sequence:
            self.logger.warning(f'ID "{self.sequence_id}" failed on BsaI restriction sites')
        else:
                self.logger.info("pass")
        self.logger.info(f'Checking BsmBI restriction sites')
        if 'CGTCTC' in self.sequence or 'GAGACG' in self.sequence:
            self.logger.warning(f'ID "{self.sequence_id}" failed on BsmBI restriction sites')
        else:
                self.logger.info("pass")
        self.logger.info(f'Checking too short sequence')
        if len(self) < 10:
            self.logger.warning(f'ID "{self.sequence_id}" failed on small sequence (<10 bp)')
        else:
                self.logger.info("pass")
        self.logger.info(f'Checking too long sequence')
        if len(self) > 20000:
            self.logger.warning(f'ID "{self.sequence_id}" failed on long sequence (>20 kbp)')
        else:
                self.logger.info("pass")
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
        self.logger.warning(f'{codon} --> {new_codon}')
        if '|REC' not in self.sequence_id:
            self.sequence_id += '|REC'
        self.sequence = f'{self.sequence[:i]}{new_codon}{self.sequence[i+3:]}'

        return self

    def remove_cutsites(self, restriction_enzymes, table=DEFAULT_TABLE):
        '''Remove recognition sites for restriction enzymes.'''
        self.logger.info(f'Removing cutsites for {restriction_enzymes}')
        changes = 0
        for enzyme in restriction_enzymes:
            for cutsite in enzyme.cutsite_list:
                while cutsite in self.sequence:
                    self.logger.warning(f'{enzyme.enzyme_id} cuts ({cutsite})')
                    changes += 1
                    self = self.recode_sequence(cutsite, table=table)
        self.logger.info(f'remove_cutsites made {changes} changes in {self.sequence_id} sequence')

        return self

    def optimize_codon_usage(self, table=DEFAULT_TABLE, maximum=False):
        '''Optimize codon usage of a given DNA sequence'''
        self.logger.debug('Optimizing codon usage...')
        if not self.basic_cds:
            self = bool_user_prompt(self, 'Optimize')
            if 'FORCED' not in self.sequence_id:
                return self
        seq_id = self.sequence_id
        optimized = self.translate(table=table).reverse_translate(table=table, maximum=maximum)

        return Nucleotide(f'{seq_id}|OPT', optimized.sequence)

    def make_part(self, part_type='3t', part_options=GGA_PART_TYPES):
        '''Make DNA part out of a given sequence'''
        self.logger.debug('Making parts...')
        seq_id = f'part_gge{part_type}_{self.sequence_id}'
        part = part_options[f'type{part_type}']
        if part_type in ('3t', '3a', '3b') and self.translate(check=True).sequence[-1] == '*':
            sequence = f'{part["prefix"]}{self.sequence[:-3]}{part["suffix"]}'
        else:
            sequence = f'{part["prefix"]}{self.sequence}{part["suffix"]}'

        return Nucleotide(seq_id, sequence)
    
    def harmonize(self, source_table, mode=0, table=DEFAULT_TABLE):
        '''Optimize codon usage of a given DNA sequence
        mode: 0 for closest frequency; 1 for same index'''
        self.logger.debug('Starting special optimization using source codon table...')
        if not self.basic_cds:
            self = bool_user_prompt(self, 'Special optimize')
            if 'FORCED' not in self.sequence_id:
                return self
        
        self.logger.info(f'Special optimization: {"closest frequency" if mode == 0 else "same index"}')

        seq_id = self.sequence_id
        optimized = list()

        for amino, triplet in zip(self.translate(table=table).sequence, self.make_triplets()):
            if amino == '?':
                self.logger.warning(f'Unknown amino acid! Appending "NNN"')
                optimized.append('NNN')
            else:
                codons = table.loc[amino]
                source_codons = source_table.loc[amino]
                sorted_codons_freq = sorted(codons['Fraction'])
                source_codon_freq = source_codons.loc[triplet]['Fraction']

                if mode == 0:
                    best, freq = 1, 0
                    for cod in sorted_codons_freq:
                        current_best = abs(cod - source_codon_freq)
                        if current_best < best:
                            best, freq = current_best, cod

                    closest_freq_codon = codons[codons['Fraction'] == freq].index[0]
                    self.logger.info(f'{triplet} (f:{source_codon_freq:.3f}) --> {closest_freq_codon} (f:{freq:.3f})')
                    optimized.append(closest_freq_codon)
                
                elif mode == 1:
                    sorted_source_codons = sorted(source_codons['Fraction'])
                    source_codon_index = sorted_source_codons.index(source_codons.loc[amino]['Fraction'])
                    same_index_codon = codons[codons['Fraction'] == sorted_codons_freq[source_codon_index]].index[0]
                    self.logger.info(f'{triplet} --> {same_index_codon} (i:{source_codon_index})')
                    optimized.append(same_index_codon)
                
                else:
                    self.logger.error('Unsupported mode, expected: 1 or 0. Skipping operation')
                    return self
        
        return Nucleotide(f'{seq_id}|SOPT{mode}', ''.join(optimized))
    
    def data_fraction(self, table=DEFAULT_TABLE, window=16):
        '''Calculates average window codon fraction for a given sequence and codon usage table.
        Returns a list of window-fraction values, which can be used for analysis or ploted.'''

        values, data = [], []
        codons = table.reset_index().set_index(['Triplet'])

        for triplet in self.make_triplets():
            values.append(codons.loc[triplet]['Fraction'])

        for n in range(len(values)+1-window):
            data.append(sum([f for f in values[n:n+window]]) / window)

        return data
    
    def data_minmax(self, table=DEFAULT_TABLE, window=16):
        '''Calculates the %MinMax values for a given sequence and codon usage table.
        Returns a list of %MinMax values, which can be used for analysis or ploted.

        Reference:
        Clarke TF IV, Clark PL (2008) Rare Codons Cluster. PLoS ONE 3(10): e3412.
        doi:10.1371/journal.pone.0003412'''

        tri_table = table.reset_index(level='Triplet')
        values, data = [], []

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

            if other_id:
                target = species

            ax1.plot(x, y2, alpha=0.8, linewidth=.5)
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
        plt.xlabel('Codon No.')
        plt.show()

        return


class Enzyme:
    def __init__(self, enzyme_id, substrate, enzyme_description=None, logger=None):
        self.enzyme_id = enzyme_id
        self.substrate = substrate
        self.enzyme_description = enzyme_description
        self.logger = logger or logging.getLogger(__name__)

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
        self.logger.debug('Generating cutsite_list from Restriction_Enzyme')
        return [self.recognition_sequence.sequence, self.recognition_sequence.reverse_complement.sequence]
