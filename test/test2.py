import pandas as pd
from os import path
from pprint import pprint

codon_table = pd.read_csv("../seqtools/data/sample_table.csv", header=None)

cds = "MDYNSADFKEIWGKAADTALLGPYNYLANNRGHNIREHLIAAFGAVIKVDKSDLETISHITKILHNSSLLVDDVEDNSMLRRGLPAAHCLFGVPQTINSANYMYFVALQEVLKLKSYDAVSIFTEEMINLHRGQGMDLYWRETLTCPSEDEYLEMVVHKTGGLFRLALRLMLSVASKQEDHEKINFDLTHLTDTLGVIYQILDDYLNLQSTELTENKGFCEDISEGKFSFPLIHSIRTNPDNHEILNILKQRTSDASLKKYAVDYMRTETKSFDYCLKRIQAMSLKASSYIDDLAAAGHDVSKLRAILHYFVSTSDCEERKYFEDAQ*"

def reverse_translate(amino, codon_table):
    """
    Reverse-translates aminoacid back to the most used codon from codon table.
    """
    temp_values = []
    temp_triplets = []

    for i, row in codon_table.loc[codon_table[0] == amino].iterrows():
        temp_triplets.append((row[1]).strip().upper())
        temp_values.append(row[2])
    
    return temp_triplets[temp_values.index(max(temp_values))]


def load_codon_usage_table(codon_table):
    # print(codon_table)

    codon_usage_aa = {}
    codon_usage_nuc = []

    for x, y in codon_table.iterrows():
        if y[0] not in codon_usage_aa:
            codon_usage_aa[y[0]] = {}
        codon_usage_aa[y[0]][y[1]] = [float(y[2])]


        codon_usage_nuc.append([y[1], y[0], y[2]])

    pprint(codon_usage_aa)
    print()
    pprint(sorted(codon_usage_nuc, key=lambda x: x[1]))




load_codon_usage_table(codon_table)



# print("".join([reverse_translate(a, codon_table) for a in cds]))