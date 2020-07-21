""""
Считаем в базе количество повторяющихся участков в последовательности
"""
from collections import Counter
from itertools import groupby
from Bio import SeqIO

AA_LENGTH = 4 # минимальная длина последовательности

final_lst = []

count = set() # id fasta

def long_repeat(line):
    """
        return the substrings that consists of the same char
    """
    lst_substr = []
    gr = [list(g) for k,g in groupby(line)]
    return list(filter(lambda x: len(x)>=AA_LENGTH, gr))
      

# DB file
file_input = r"C:\Users\Фесенко\Google Диск\lncRNAs_sORFs\lncRNAs_sorfs_mipepid\combined_sORFs_mipepid_locus_transcripts_nonredundant_nonested.fasta"


#reading DB
for record in SeqIO.parse(file_input, 'fasta'):
    res = long_repeat(str(record.seq))
    if len(res)>0:
        count.add(record.id)
    final_lst.extend(["".join(i) for i in res])
    
# out data
with open(r"C:\Users\Фесенко\Google Диск\lncRNAs_sORFs\Paper\aa_string_count4.txt", 'w') as out:
    for i in final_lst:
        out.write(f"{i}\n")

#Counting data
c = Counter(final_lst)

print(c.most_common(10))
print(len(count))