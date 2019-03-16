'''Tests for modules Sequence, Protein and Nucleotide.'''
import pytest
import pandas as pd
from seqtools.modules import Sequence, Protein, Nucleotide


def bio_sequence(name, seq):
    return Sequence(name, seq)


def protein_sequence(name, seq):
    return Protein(name, seq)


def nucleotide_sequence(name, seq):
    return Nucleotide(name, seq)


def test_bio():
    bio = bio_sequence('bio_01', 'abcdefghabcdefghabcdefghabcdefgh')
    assert bio.sequence_id == 'bio_01'
    assert bio.sequence != 'abcdefghabcdefghabcdefghabcdefgh'
    assert bio.sequence == 'abcdefghabcdefghabcdefghabcdefgh'.upper()
    assert len(bio) == len('abcdefghabcdefghabcdefghabcdefgh')
    assert str(bio)[:7] == '>bio_01'


def test_wrong_protein():
    with pytest.raises(ValueError):
        protein_sequence('prot_00', 'abcdefgabcdefgabcdefgabcdefg')

    with pytest.raises(ValueError):
        protein_sequence('prot_00', '123')


def test_wrong_nucleotide():
    with pytest.raises(ValueError):
        nucleotide_sequence('dna_00', 'actgeactgeactgeactgeactgeactgeactge')


def test_correct_protein():
    prot = protein_sequence('prot_01', 'GALMFWKQESPVICYHRNDTGALMFWKQESPVICY?*')
    assert prot.sequence_id == 'prot_01'
    assert prot.sequence != 'GALMFWKQESPVICYHRNDTGALMFWKQESPVICY?*'.lower()
    assert prot.sequence == 'GALMFWKQESPVICYHRNDTGALMFWKQESPVICY?*'
    assert len(prot) == len('GALMFWKQESPVICYHRNDTGALMFWKQESPVICY?*')
    assert str(prot)[:8] == '>prot_01'
    new = prot + prot
    assert new.sequence_id == 'concatenated'
    assert new.sequence == 'GALMFWKQESPVICYHRNDTGALMFWKQESPVICY?*GALMFWKQESPVICYHRNDTGALMFWKQESPVICY?*'
    assert len(new) == len('GALMFWKQESPVICYHRNDTGALMFWKQESPVICY?*GALMFWKQESPVICYHRNDTGALMFWKQESPVICY?*')


def test_correct_nucleotide():
    dna = nucleotide_sequence('dna_01', 'actgactgactgactgactgactgn')
    rna = nucleotide_sequence('rna_01', 'acugacugacugacugacugacugn')
    assert dna.sequence_id == 'dna_01'
    assert rna.sequence_id == 'rna_01'
    assert dna.sequence != 'actgactgactgactgactgactgn'
    assert dna.sequence == 'actgactgactgactgactgactgn'.upper()
    assert len(dna) == len('actgactgactgactgactgactgn')
    assert str(dna)[:7] == '>dna_01'
    assert dna.is_cds is False
    new = dna + dna
    assert new.sequence_id == 'concatenated'
    assert new.sequence == 'ACTGACTGACTGACTGACTGACTGNACTGACTGACTGACTGACTGACTGN'
    assert len(new) == len('ACTGACTGACTGACTGACTGACTGNACTGACTGACTGACTGACTGACTGN')


dummy_table = pd.read_csv('data/sample_table.csv', header=None)


def test_protein_operations():
    prot = protein_sequence('prot_02', 'GALMFWKQESPVICYHRNDTGALMFWKQESPVICY*')
    assert prot.sequence_id == 'prot_02'
    assert prot.sequence == 'GALMFWKQESPVICYHRNDTGALMFWKQESPVICY*'
    dna = prot.reverse_translate(dummy_table)
    assert dna.sequence_id == 'prot_02|NUC'
    assert dna.sequence == 'GGCGCCCTCATGTTCTGGAAGCAGGAGTCCCCAGTGATCTGCTACCACCGAAACGACACCGGCGCCCTCATGTTCTGGAAGCAGGAGTCCCCAGTGATCTGCTACTAA'


def test_dna_operations():
    dna = nucleotide_sequence('dna_02', 'ATGGCTCGAAAGGCTCCTCACATCGACTAA')
    assert dna.sequence_id == 'dna_02'
    assert dna.sequence == 'ATGGCTCGAAAGGCTCCTCACATCGACTAA'
    prot = dna.translate(dummy_table)
    assert prot.sequence_id == 'dna_02|PROT'
    assert prot.sequence == 'MARKAPHID*'
    opt = dna.optimize_codon_usage(dummy_table)
    assert opt.sequence_id == 'dna_02|OPT'
    assert opt.sequence == 'ATGGCCCGAAAGGCCCCACACATCGACTAA'
