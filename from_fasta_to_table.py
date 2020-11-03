# получаем сиквенсы и имена из fasta в таблицу

from Bio import SeqIO
import pandas as pd

# read file
file_input = r"/Users/igorfesenko/Google Диск/lncRNAs_sORFs/DB/combined_sORFs_mipepid_locus_transcripts_nonredundant_nonested90317.fasta"

# out table
out_file = r"/Users/igorfesenko/Google Диск/lncRNAs_sORFs/DB/MiPepid90317_id_seq.xlsx"

results = []

for record in SeqIO.parse(file_input, 'fasta'):
    results.append([record.id, record.seq, len(record.seq), record.description])

df = pd.DataFrame(columns=['ID', 'Sequence', 'Length', 'Description'], data=results)

df.to_excel(out_file, index=False)