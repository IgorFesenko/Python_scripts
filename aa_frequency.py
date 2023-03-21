"""
Считаем частоту каждой аминокислоты в базе данных последовательностей
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import Counter
import pandas as pd

#db_file = r'/Users/igorfesenko/Google Диск/lncRNAs_sORFs/DB/filt_proteincoding_proteins.fasta_asterix.fasta'
#db_file = r'/Users/igorfesenko/Google Диск/lncRNAs_sORFs/lncRNAs_sorfs_mipepid/combined_sORFs_mipepid_locus_transcripts_nonredundant_nonested90317.fasta'
#db_file = r'/Users/igorfesenko/filt_proteincoding_proteins.fasta_asterix_small.fasta'
#db_file = r'/Users/igorfesenko/Google Диск/lncRNAs_sORFs/DB/filt_proteincoding_proteins16178_asterix.fasta'
#db_file = r'/Users/igorfesenko/Google Диск/lncRNAs_sORFs/DB/filt_proteincoding_proteins_small252asterix.fasta'
db_file = r'/Users/igorfesenko/Google Диск/lncRNAs_sORFs/DB/ilt_proteincoding_proteins_big15926asterix.fasta'



aa_lst = []


for record in SeqIO.parse(db_file, 'fasta'):
    lst = list(str(record.seq))
    aa_lst.extend(lst)

print(f"The number of AA in db: {len(aa_lst)}")

c = Counter(aa_lst)

print(c)

print(c.most_common(10))

df = pd.DataFrame(columns=['aa','number'], data=c.most_common(25))

df['frequency'] = df['number'].map(lambda x: x/len(aa_lst))

print(df)
df.to_csv(r'/Users/igorfesenko/Google Диск/lncRNAs_sORFs/LCR/big_proteins15926_aa_freq.csv', index=False)
