"""
Извлекаем длину мРНК, которые используем для сравнения
"""

from Bio import SeqIO
import pandas as pd

fasta_input = r'/Users/igorfesenko/Google Диск/Phytozome_Ppatens_29012020/PhytozomeV12/Ppatens/annotation/Ppatens_318_v3.3.transcript_primaryTranscriptOnly.fa'

defline = r'/Users/igorfesenko/Google Диск/Phytozome_Ppatens_29012020/PhytozomeV12/Ppatens/annotation/Ppatens_318_v3.3.defline.txt'

out_file = r'/Users/igorfesenko/Google Диск/lncRNAs_sORFs/Paper/mRNAs_length.txt'

df = pd.read_csv(defline, sep='\t', names=['prot','pdef','desc'])

prot_lst = set(df['prot'])

print(len(prot_lst))

# indexing of sequences
dict_db1 = SeqIO.index(fasta_input, "fasta")

with open (out_file, 'a') as out:
    for p in prot_lst:
        out.write(f"{p} {len(dict_db1[p].seq)}\n")

    




#print(df.head())