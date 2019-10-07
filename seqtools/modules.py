#cSpell: Disable#
import re
from seqtools.util import bool_user_prompt, sequence_match, get_codon


class Sequence:
    '''Biological sequence object'''

    def __init__(self, sequence_id, sequence):
        self.sequence_id = sequence_id
        self.sequence = sequence.upper()

    def __repr__(self):
        return f'{self.sequence_id}, {self.sequence}'

    def __str__(self):
        return f'>{self.sequence_id}\n{self.sequence}'

    def __len__(self):
        return len(self.sequence)
    
    ### check how to check if seq is prot or nuc
    ### then you can delete __add__ methods from the Protein and Nucleotide classes
    # def __add__(self, other):
    #     return self.sequence + other.sequence

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

    def reverse_translate(self, table, maximum=False):
        '''Returns optimized DNA sequence'''
        dna_sequence = list()
        if maximum:
            name = f'{self.sequence_id}|NUC-MAX'
            for amino in self.sequence:
                dna_sequence.append(get_codon(table, amino, maximum=True))
        else:
            name = f'{self.sequence_id}|NUC'
            for amino in self.sequence:
                dna_sequence.append(get_codon(table, amino))

        return Nucleotide(name, ''.join(dna_sequence))


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
    def is_cds(self):
        '''Returns True if sequence is CDS or false if its not'''
        if self.sequence[:3] == 'ATG' and len(self) % 3 == 0:
            return True
        else:
            return False

    def reverse_complement(self):
        '''Returns reverse complement of given DNA sequence'''
        return Nucleotide(f'{self.sequence_id}|REVCOMP',
            self.sequence.translate(str.maketrans("ACGT", "TGCA")[::-1]))

    def make_triplets(self):
        '''Makes list of chunks 3 characters long from a sequence'''
        return [self.sequence[start:start + 3] for start in range(0, len(self.sequence), 3)]

    def translate(self, table):
        '''Translate DNA sequence in PROTEIN sequence'''
        prompt = f'Sequence with ID {self.sequence_id} is not a CDS. Translate anyway?'

        if not self.is_cds:
            if not bool_user_prompt(prompt):
                return None
            else:
                new_sequence_id = f'{self.sequence_id}|FORCED'
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
    
    # def generate_codon_dictionary(self, codon_table):
    #     '''Returns dictionary to store the codons that have been attempted thus far'''
    #     return dict(enumerate(zip(self.translate(codon_table),
    #                 [[self.sequence[n:n+3]] for n in range(0, self.len, 3)])))

    def optimize_codon_usage(self, codon_table):
        '''Optimize codon usage given with codon usage table'''
        prompt = f'Sequence with ID {self.sequence_id} is not a CDS. Optimize anyway?'

        if not self.is_cds:
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
        '''
        if not self.is_cds:
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

        return sum(values) / len(values)
