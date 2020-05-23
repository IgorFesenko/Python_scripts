
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re

db_file = r"C:\Users\fesenkoi2\IFESENKO\lncRNAs_sORFs\lncRNAs_sorfs_mipepid\combined_sORFs_mipepid_locus_transcripts_nonredundant.fasta"

dict_db = {}

for record in SeqIO.parse(db_file, 'fasta'):
    dict_db[str(record.seq)]=record.id


cnt=0    
for record in SeqIO.parse(db_file, 'fasta'):
    search_space = {k:v for k,v in dict_db.items() if len(k)>len(record.seq)}
    #print(search_space.keys())
    for k,v in search_space.items():
        if re.findall(str(record.seq), k):
            cnt+=1
            print(cnt)
            with open (r"C:\Users\fesenkoi2\IFESENKO\lncRNAs_sORFs\lncRNAs_sorfs_mipepid\nested_combined_sORFs_mipepid.txt", 'a') as out_names:
                out_names.write("{} {}".format(record.id,v)+'\n')
            #print(record.seq, 'in', k)

print(f"The number of nested sORFs: {cnt}")