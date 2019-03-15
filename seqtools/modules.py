from seqtools.util import make_triplets, bool_user_prompt


class Sequence:
    '''
    Biological sequence object.
    '''

    def __init__(self, sequence_id, sequence):
        self.sequence_id = sequence_id
        self.sequence = sequence.upper()

    def __repr__(self):
        return f'{self.sequence_id}, {self.sequence}'

    def __str__(self):
        return f'>{self.sequence_id}\n{self.sequence}'

    def __len__(self):
        return len(self.sequence)

    @property
    def triplets(self):
        return make_triplets(self.sequence)


class Protein(Sequence):
    '''
    PROTEIN sequence object.
    '''
    allowed_characters = '*?galmfwkqespvicyhrndt'

    def __init__(self, sequence_id, sequence):
        super().__init__(sequence_id, sequence)

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
    '''
    DNA/RNA sequence object.
    '''

    allowed_characters = 'atugcn'

    def __init__(self, sequence_id, sequence):
        super().__init__(sequence_id, sequence)

    def translate(self, codon_table):
        '''
        Translate DNA in PROTEIN.
        '''
        if self.sequence[:3] != 'ATG' or len(self) % 3 != 0:
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
        '''
        Optimize codon usage given with codon usage table.
        '''
        if len(self) % 3 != 0:
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
