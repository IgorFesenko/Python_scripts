# фильтрация fasta по заданному списку

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re

# file with initial sequences 
file_input = r"C:\Users\fesenkoi2\IFESENKO\orthologs\Ppatens_318_v3.3.cds.fa"

# out file
records = []
out_file = r"C:\Users\fesenkoi2\IFESENKO\orthologs\Ppatens_318_v3.3.cds_filt.fa"

# list of names
file_filt = r'C:\Users\fesenkoi2\IFESENKO\filt_proteincoding_proteins_id.txt'

#reading ID from file
with open(file_filt) as inp:
    lst = set([i.strip()[:-2] for i in inp.readlines()])


print(f"The number of records: {len(lst)}")

# indexing of sequences
dict_db1 = SeqIO.index(file_input, "fasta")

# creating new fasta file
for sorf in lst:
    records.append(SeqRecord(seq=dict_db1[sorf].seq, id=sorf))

SeqIO.write(records,out_file,"fasta")