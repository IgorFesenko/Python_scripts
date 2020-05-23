#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Running PAL2NAL
'''
import os
import sys
import getopt
import pandas as pd
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO



if __name__ == "__main__":
   
    unixOptions = "f:i:"  
    gnuOptions = ["file=","id="]
    # Получаем строку входных параметров
    fullCmdArguments = sys.argv
    argumentList = fullCmdArguments[1:]

    try:  
        arguments, values = getopt.getopt(argumentList, unixOptions, gnuOptions)
    except getopt.error as err:  
        print (str(err))
        sys.exit(2)      

    
    file = '' 
    ID = '' 
    for currentArgument, currentValue in arguments:  
        if currentArgument in ("-f", "--file"):
            file = currentValue
        elif currentArgument in ("-i", "--id"):
            ID = currentValue                                   
    

    prot_file = r"./protein.fasta"
    nucl_file = r"./nucleotide.fasta"


    df = pd.read_csv(file, compression='gzip')

    # fitration by organism
    res_table = df[df['id']==ID]

    #remove duplucates by aligment
    res_table = res_table.sort_values(by='align_len', ascending=False)
    res_table.drop_duplicates(subset=['sORF'], keep='first', inplace=True)

    for index,row in res_table.iterrows():
        query_pep = Seq(row['Seq_query'], IUPAC.protein)
        hit_pep = Seq(row['Seq_hit'], IUPAC.protein)

        prot=[SeqRecord(seq=query_pep, id=row['sORF']),SeqRecord(seq=hit_pep, id=row['hit_id'])]
        
        query_nucl = Seq(row['query_nucleotide'], IUPAC.unambiguous_dna)
        hit_nucl = Seq(row['hit_nucleotide'], IUPAC.unambiguous_dna)
        
        nucl = [SeqRecord(seq=query_nucl, id=row['sORF']), SeqRecord(seq=hit_nucl, id=row['hit_id'])]

        SeqIO.write(prot,prot_file,"fasta")
        SeqIO.write(nucl,nucl_file,"fasta")

        os.system("perl ./pal2nal.pl {} {} -output fasta > ./tmp/file{}.fasta".format(prot_file, nucl_file,index))

    