#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Running MAFFT and save alignments
'''

import os
import sys
import getopt
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

if __name__ == "__main__":

    # reading files
    cwd = os.getcwd()
    names = [n for n in os.listdir() if n.endswith('.fasta')]



    # creating directory
    try:
        os.mkdir('mafft_align')
    except:
        os.system('rm -r ./mafft_align')
        os.mkdir('mafft_align')

    # running mafft
    for n in names:
        print(n)
        os.system(f"/usr/local/bin/mafft --thread 4 --localpair  --maxiterate 16 --phylipout --inputorder {n}  > ./mafft_align/{n.split('.')[0]}.phy")

