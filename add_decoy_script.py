
"""
Добавляем decoy
"""

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

names = [i for i in os.listdir() if i.endswith('fasta')]

for name in names:
    records = []
    records_rev = []
    # indexing of sequences
    #dict_db = SeqIO.index(name, "fasta")
    for record in SeqIO.parse(name, 'fasta'):
        records.append(SeqRecord(seq=record.seq, id=record.id, description=record.description))
        records_rev.append(SeqRecord(seq=record.seq[::-1], id=f"rev_{record.id}"))
    SeqIO.write(records+records_rev,f"{name}_decoy.fasta","fasta")
