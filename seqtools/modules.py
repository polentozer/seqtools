import re
from seqtools.util import make_triplets, bool_user_prompt, sequence_match


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

    def reverse_translate(self, codon_table):
        reverse_translation = ''
        for amino in self.sequence:
            temp_values = []
            temp_triplets = []

            for _, row in codon_table.loc[codon_table[0] == amino].iterrows():
                temp_triplets.append(row[1])
                temp_values.append(row[2])

            reverse_translation += temp_triplets[temp_values.index(max(temp_values))]

        return Nucleotide(f'{self.sequence_id}|NUC', reverse_translation)


class Nucleotide(Sequence):
    '''DNA/RNA sequence object'''

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
    def triplets(self):
        return make_triplets(self.sequence)

    @property
    def is_cds(self):
        '''Returns True if sequence is CDS or false if its not'''
        if self.sequence[:3] == 'ATG' and len(self) % 3 == 0:
            return True
        else:
            return False

    def translate(self, codon_table):
        '''Translate DNA sequence in PROTEIN sequence'''

        if not self.is_cds:
            if not bool_user_prompt(f'Sequence with ID {self.sequence_id} is not a CDS. Translate anyway?'):
                return None
            else:
                new_sequence_id = f'{self.sequence_id}|FORCED'
        else:
            new_sequence_id = self.sequence_id

        translation = ''
        for triplet in self.triplets:
            if len(triplet) == 3:
                translation += codon_table.loc[codon_table[1] == triplet][0].iloc[0][0]
            else:
                translation += '?'

        return Protein(f'{new_sequence_id}|PROT', translation)

    def optimize_codon_usage(self, codon_table):
        '''Optimize codon usage given with codon usage table'''

        if not self.is_cds:
            if not bool_user_prompt(f'Sequence with ID {self.sequence_id} is not a CDS. Optimize anyway?'):
                return self
            else:
                new_sequence_id = f'{self.sequence_id}|FORCED'
        else:
            new_sequence_id = self.sequence_id

        optimized_sequence = ''
        for triplet in self.triplets:
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
        pass
