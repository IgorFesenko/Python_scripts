"""
Формируем fasta файл из таблицы полученной при работе скрипта full_seq_extractor
"""
import pandas as pd
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np



table = pd.read_excel(r"F:\YandexDisk-feigor\1kb_sorfs_logo\ORF19_TDM\script_results\table_seq.xlsx")
filt_table = table[np.logical_not(table['Full_PEP']==None)]
out_fasta = r"F:\YandexDisk-feigor\1kb_sorfs_logo\ORF19_TDM\script_results\full_seq.fasta"

results=[]

for index,row in filt_table.iterrows():
    seq_id = f"{row['ID']}_{row['species']}"
    if row['Full_PEP'] != 'None':
        #print(row['Full_PEP'])
        seq_seq = Seq(row['Full_PEP'], IUPAC.unambiguous_dna)
        results.append(SeqRecord(seq=seq_seq, id=seq_id))
    else:
        continue
        

       

SeqIO.write(results,out_fasta,"fasta")