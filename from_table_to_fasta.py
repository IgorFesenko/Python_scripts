"""
Формируем fasta файл из таблицы полученной при работе скрипта full_seq_extractor
"""
import pandas as pd
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO



table = pd.read_excel(r"F:\YandexDisk-feigor\1kb_sorfs_logo\DVL_peptides\Pp3c19_10140V3.1_ORF20\table_seq.xlsx")
out_fasta = r"F:\YandexDisk-feigor\1kb_sorfs_logo\DVL_peptides\Pp3c19_10140V3.1_ORF20\full_seq.fasta"

results=[]

for index,row in table.iterrows():
    seq_id = row['ID']
    if row['Full_PEP'] != None:
        seq_seq = Seq(row['Full_PEP'], IUPAC.unambiguous_dna)
        results.append(SeqRecord(seq=seq_seq, id=seq_id))
    else:
        continue
        

       

SeqIO.write(results,out_fasta,"fasta")