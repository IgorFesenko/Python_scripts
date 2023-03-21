#!/Users/igorfesenko/anaconda3/bin/python3
# -*- coding: utf-8 -*-

"""
Creating fasta files with nucleotides and proteins smORFs

These files will be analyzed by MAFFT, IQ-TREE and HyPHY
"""

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import getopt
import pandas as pd
from Bio.Seq import Seq





if __name__ == "__main__":
   
    #reading arguments    
    unixOptions = "f:"  
    gnuOptions = ["file="]
    #arguments
    fullCmdArguments = sys.argv
    argumentList = fullCmdArguments[1:]

    try:  
        arguments, values = getopt.getopt(argumentList, unixOptions, gnuOptions)
    except getopt.error as err:  
        print (str(err))
        sys.exit(2)      

    
    file = '' 
   
    for currentArgument, currentValue in arguments:  
        if currentArgument in ("-f", "--file"):
            file = currentValue


    # creating directory
    try:
        os.mkdir('tmp_fasta')
    except:
        os.system('rm -r tmp_fasta')
        os.mkdir('tmp_fasta')
    

    #reading table
    df = pd.read_csv(file, compression='gzip')

    # reading column names
    prot_names = [i for i in df.columns if 'aa' in i]
    nucl_names = [i for i in df.columns if 'nucl' in i]
    id_names = [i for i in df.columns if 'hit' in i]

    # Parsing table

    for index,row in df.iterrows():
        file_name = row['Protein']
        prot_results = []
        nucl_results = []
        # protein fasta
        prot_results.append(SeqRecord(seq=Seq(row['Sequence']), id=row['sORF'], description=""))
        nucl_results.append(SeqRecord(seq=Seq(row['Sequence_n']), id=row['sORF'], description=""))
        for n in prot_names:
            if row[n]:
                hit_name = f"hit_{n.split('_')[1]}"
                prot_results.append(SeqRecord(seq=Seq(row[n]), id=hit_name, description=""))
        for n in nucl_names:
            if row[n]:
                hit_name = f"hit_{n.split('_')[1]}"
                nucl_results.append(SeqRecord(seq=Seq(row[n]), id=hit_name, description=""))
    
        SeqIO.write(prot_results,f"./tmp_fasta/PROT_{file_name}.fasta","fasta")
        SeqIO.write(nucl_results,f"./tmp_fasta/PROT_{file_name}.fasta","fasta")
            


