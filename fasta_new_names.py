from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd


# creating fasta with new names from paper

# read file
file_input = r'/Users/igorfesenko/Google Диск/lncRNAs_sORFs/DB/smORFs_mipepid90057.fasta'

# out file
records = []# new SeqRecords
out_file = r'/Users/igorfesenko/Google Диск/lncRNAs_sORFs/DB/smORFs_mipepid90057_paper_name.fasta'

# table names 

names_table = pd.read_excel(r'/Users/igorfesenko/Google Диск/lncRNAs_sORFs/Paper/suppl_tables/mipepid_table_allsorfs_filt90057name.xlsx')



# parsing data
for record in SeqIO.parse(file_input, 'fasta'):
    new_name = names_table[names_table['sORF']==str(record.id)]['name'].values[0]
    #print(new_name)
    records.append(SeqRecord(seq=record.seq, id=str(new_name), description=record.description))
    
    

SeqIO.write(records,out_file,"fasta")
