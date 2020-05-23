#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Creating two fasta from PAL2NAL output
'''
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

if __name__ == "__main__":
    names = [n for n in os.listdir() if n.endswith('.fasta')]

    records_physco = []
    records_other = []

    for name in names:
        tmp_fasta =[]
        for record in SeqIO.parse(name, 'fasta'):
            tmp_fasta.append(SeqRecord(seq=record.seq, id=record.id, description=record.description))

        records_physco.append(tmp_fasta[0])
        records_other.append(tmp_fasta[1])

    
    SeqIO.write(records_physco,'./physco_orthologs.fasta',"fasta")

    SeqIO.write(records_other,'./other_species.fasta',"fasta")