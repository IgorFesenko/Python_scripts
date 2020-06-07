#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Running psi-blast
"""

import os
import sys
import getopt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


if __name__ == "__main__":
    
    unixOptions = "q:d:f"  
    gnuOptions = ["query=", "db=","outfmt="]
    
    fullCmdArguments = sys.argv
    argumentList = fullCmdArguments[1:]

    
    try:  
        arguments, values = getopt.getopt(argumentList, unixOptions, gnuOptions)
    except getopt.error as err:  
        print (str(err))
        sys.exit(2)      

    query = ''
    db = ''
    outfmt = ''   
    for currentArgument, currentValue in arguments:  
        if currentArgument in ("-q", "--query"):
            query = currentValue                                   
        elif currentArgument in ("-d", "--db"):
            db = currentValue
        elif currentArgument in ("-f", "--outfmt"):
            outfmt = currentValue

    print(query,db,outfmt)       
    os.system("mkdir tmp_psiblast") # creating tmp directory

    #creating fragments of query file
    names = set()
    def fasta_record(rec,cnt):
        name = './tmp_psiblast/fragment_{}.fasta'.format(cnt)
        names.add(name)
        print(name)
        SeqIO.write(rec,name,"fasta")


    seq_count=0
    records=[]
    cnt=0
    for record in SeqIO.parse(query, 'fasta'):
        cnt+=1
        seq_count+=1
        if seq_count<100:
            records.append(SeqRecord(seq=record.seq, id=record.id, description=record.description))
            
        else:
            fasta_record(records,cnt)
            records=[]
            seq_count=0
    fasta_record(records,cnt)
    print('The number of records: ',cnt)

    # running blast
    for n in names:
        print('running..',n)
        out_name = "{}_psiblastx_{}".format(n , db.split('/')[-1])
        os.system("/usr/bin/psiblast -db {} -query {} -outfmt {} -out {} &".format(db,n,outfmt,out_name))