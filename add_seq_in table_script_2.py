#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import getopt
import pandas as pd
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from time import process_time
import datetime
import swifter


if __name__ == "__main__":
    unixOptions = "t:q:h:"  
    gnuOptions = ["table=","db1=","db2="]

    fullCmdArguments = sys.argv
    argumentList = fullCmdArguments[1:]

    try:  
        arguments, values = getopt.getopt(argumentList, unixOptions, gnuOptions)
    except getopt.error as err:  
        print (str(err))
        sys.exit(2)      

    table = ''
    db1 = ''
    db2 = ''
    
       
    for currentArgument, currentValue in arguments:  
        if currentArgument in ("-t", "--table"):
            table = currentValue
        elif currentArgument in ("-q", "--db1"):
            db1 = currentValue
        elif currentArgument in ("-h", "--db2"):
            db2 = currentValue 


    b = datetime.datetime.now()
    print(f"The program started at {b.strftime('%Y-%m-%d %H.%M')}")


    def query_sequence(row):
        start = int(row['query_start'])*3-3
        stop = int(row['query_stop'])*3
        seq_sorf = str(dict_db2[row['Protein']].seq)
        return seq_sorf[start:stop]


    def hit_sequence(row):
        start = int(row['hit_start'])-1
        stop = int(row['hit_stop'])
        seq_transcipt = str(dict_db1[str(row['hit_id'])].seq)
        seq_seq = Seq(seq_transcipt[start:stop], IUPAC.unambiguous_dna)
        if '-'in row['frame']:
            return [row['hit_id'],str(seq_seq.reverse_complement())]
        else:
            return [row['hit_id'],seq_seq]

    print('READING TABLE...')

    table_file=table

    df = pd.read_csv(table_file, compression='gzip')


    print('READ DATABASE...')
    t_start2 = process_time()

    dict_db1 = SeqIO.index(db1, "fasta") # mosses transcripts
    dict_db2 = SeqIO.index(db2, "fasta") # mipepid


    print("ADD NUCLEOTIDES...")
    print('adding queries...')
    df['query_nucleotide'] = df.swifter.apply(query_sequence,axis=1)

    print('adding hits...')
    hit_seq = df.swifter.apply(hit_sequence,axis=1)

    hit_res = pd.DataFrame(columns=['hit_id','hit_nucleotide'], data=hit_seq)
    df = df.merge(hit_res, on='hit_id', how='left')


    df.to_csv(f"{table}_nucl.csv", compression='gzip', index=False)


    t_stop2 = process_time()
    print(f"TIME: {(t_stop2-t_start2)/60} min")