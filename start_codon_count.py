#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
from collections import Counter


#input_file = r'/Users/igorfesenko/Google Диск/lncRNAs_sORFs/DB/coordinates/NCBI_finder_intersect_locus_cds.fasta'

input_file = r'/Users/igorfesenko/lncRNAs_nanopore2418_dna_nonAUG.fasta'


c = Counter()

for record in SeqIO.parse(input_file, 'fasta'):
    start_codon = str(record.seq)[:3]
    c[start_codon]+=1

print(c)