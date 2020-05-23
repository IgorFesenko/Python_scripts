'''
Creating a table with id - description pairs

'''
from Bio import SeqIO

# read file
file_input = r"C:\Users\fesenkoi2\IFESENKO\lncRNAs_sORFs\lncRNAs_sorfs_NCBIfinder\merged_lncRNAs_transcripts_filtered_nonredundant.out_30_1_plus"

out_file = r"C:\Users\fesenkoi2\IFESENKO\lncRNAs_sORFs\lncRNAs_sorfs_NCBIfinder\merged_lncRNAs_transcripts_filtered_nonredundant.out_30_1_plus_names.txt"


with open (out_file, 'a') as out:

    for record in SeqIO.parse(file_input, 'fasta'):
            
        out.write(f"{record.id} {record.description}\n")