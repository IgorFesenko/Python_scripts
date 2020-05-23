"""
Translate nucleotide sequence and create fasta with sequences
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from time import process_time
import datetime

b = datetime.datetime.now()
print(f"The program started at {b.strftime('%Y-%m-%d %H.%M')}")

t_start = process_time()


seq_file = r"C:\Users\fesenkoi2\IFESENKO\lncRNAs_sORFs\combined_set\merged_lncRNAs_transcripts_filtered.fasta"
out_file = r"C:\Users\fesenkoi2\IFESENKO\lncRNAs_sORFs\combined_set\merged_lncRNAs_transcripts_filtered)my_translation.fasta"

records = []
for record in SeqIO.parse(seq_file, 'fasta'):
    #print(record)
    messenger_rna = record.seq

    translation1 = messenger_rna.translate()
    translation2 = messenger_rna[1:].translate()
    translation3 = messenger_rna[2:].translate()
    lst_translation = [translation1,translation2,translation3]
    frame_number = ['frame1', 'frame2', 'frame3']
    for i in range(3):
        pep_numb=0
        for pep in lst_translation[i].split('*'):
            if 'M' in pep:
                position = pep.find('M')
                new_seq = pep[position:]
                if len(new_seq) >= 10:
                    pep_numb += 1
                    records.append(SeqRecord(seq=new_seq, id=f"{record.id}_{frame_number[i]}_peptide{pep_numb}", description=""))

SeqIO.write(records,out_file,"fasta")

t_stop = process_time()

print(f"TIME: {(t_stop-t_start)/60} min")