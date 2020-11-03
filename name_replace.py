from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# fasta with new names


# read file
file_input = r'/Users/igorfesenko/Google Диск/lncRNAs_sORFs/DB/sORFs_mipepid_nucleotide_db90317.fasta'

# out file
records = []
out_file = r'/Users/igorfesenko/Google Диск/lncRNAs_sORFs/DB/sORFs_mipepid_nucleotide_db90317_newname.fasta'

# new name table 
table = []
out_table = r'/Users/igorfesenko/Google Диск/lncRNAs_sORFs/DB/sORFs_mipepid_nucleotide_db90317_newnames_infoseq.txt'
cnt=1

for record in SeqIO.parse(file_input, 'fasta'):
    new_name = f"ORF_{cnt}"
    records.append(SeqRecord(seq=record.seq, id=new_name))
    table.append([record.id,new_name])
    cnt+=1

SeqIO.write(records,out_file,"fasta")

with open(out_table, 'w') as out:
    for name in table:
        out.write(f"{name[0]} {name[1]}\n")
