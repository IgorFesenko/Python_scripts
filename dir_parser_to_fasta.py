#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Парсим таблицы после работы скрипта xml_parser_web_tblastn.py
"""
import os
import sys
import getopt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC



if __name__ == "__main__":
        
    unixOptions = "d:"  
    gnuOptions = ["dir="]
   
    fullCmdArguments = sys.argv
    argumentList = fullCmdArguments[1:]

   
    try:  
        arguments, values = getopt.getopt(argumentList, unixOptions, gnuOptions)
    except getopt.error as err:  
        print (str(err))
        sys.exit(2) 

    
    directory = ''
    
    for currentArgument, currentValue in arguments:
        if currentArgument in ("-d", "--dir"):
            directory = currentValue                     

    #reading smORFs database
    db_input = r'/Users/igorfesenko/Google Диск/lncRNAs_sORFs/lncRNAs_sorfs_mipepid/combined_sORFs_mipepid_locus_transcripts_nonredundant_nonested.fasta'
    dict_db1 = SeqIO.index(db_input, "fasta")
    
    # reading directories 
    dir_names = []
    for root, dirs, files in os.walk(directory):
        for name in dirs:
            dir_names.append(os.path.join(root, name))
            
    count_hits = [] # counting number of hits for smORFs

    # creating fasta
    for d in dir_names:
        print(d)
        # reading table
        df = pd.read_excel(f"{d}/table_seq.xlsx")
        filt_table = df[np.logical_not(df['Full_PEP']==None)]
        out_fasta = f"{d}/full_seq.fasta"
        # reading smORF peptide sequence
        query = list(df['query'])[0]
        number_hits = len(set(df['ID']))
        count_hits.append([query, number_hits])

        #reading peptide sequence and name
        results=[]
        results.append(SeqRecord(seq=dict_db1[query].seq, id=query))

        #creating fasta
        
        for index,row in filt_table.iterrows():
            seq_id = f"{row['ID']}_{row['species']}"
            if row['Full_PEP'] != 'None':
                seq_seq = Seq(row['Full_PEP'], IUPAC.unambiguous_dna)
                results.append(SeqRecord(seq=seq_seq, id=seq_id))
           
                   

        SeqIO.write(results,f"{d}/fasta_alignments.fasta",'fasta')

        out_info = pd.DataFrame(columns=['query','numb_of_hits'], data=count_hits)
        out_info.to_excel(f"{directory}/info_table.xlsx", index=False)
            