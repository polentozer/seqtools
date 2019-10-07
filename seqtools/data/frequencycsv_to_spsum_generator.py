#cSpell: Disable#
import pandas

# Grab table from user
table = input("Desired table : ")

# Sets up spsum format from SPSUM_LABEL ftp://ftp.kazusa.or.jp/pub/codon/current/SPSUM_LABEL
spsum_format = "CGA CGC CGG CGT AGA AGG CTA CTC CTG CTT TTA TTG TCA TCC TCG TCT AGC AGT ACA ACC ACG ACT CCA CCC CCG CCT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT AAA AAG AAC AAT CAA CAG CAC CAT GAA GAG GAC GAT TAC TAT TGC TGT TTC TTT ATA ATC ATT ATG TGG TAA TAG TGA".split()

# Read in table as pandas dataframe
csv_data = pandas.read_csv("frequency_example.csv")
for index, row in csv_data.iterrows():
    spsum_format = [x if not x == row["codon"] else str(int(float(row[table])*100)) for x in spsum_format]

properspsum = ' '.join(spsum_format)

# Add taxid (custom_x)
taxid = input("Taxonomy id of table? : ")

# Write table
with open("custom_table.spsum", "a") as custom_table:
    custom_table.write(f'{taxid}:{table}:{str(sum(list(map(int, spsum_format))))}\n{properspsum}')

# print(f'{taxid}:{table}:{str(sum(list(map(int, spsum_format))))}\n{properspsum}')
