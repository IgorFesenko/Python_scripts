#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Parsing XML blast output
'''
import os
import sys
import getopt
from Bio.Blast import NCBIXML
import pandas as pd
from time import process_time
import datetime


if __name__ == "__main__":
    
    unixOptions = "f:e:"  
    gnuOptions = ["file=",'eval=']
   
    fullCmdArguments = sys.argv
    argumentList = fullCmdArguments[1:]

   
    try:  
        arguments, values = getopt.getopt(argumentList, unixOptions, gnuOptions)
    except getopt.error as err:  
        print (str(err))
        sys.exit(2)      
       
    for currentArgument, currentValue in arguments:  
        if currentArgument in ("-f", "--file"):
            file = currentValue                     
        if currentArgument in ("-e", "--eval"):
            evalue = float(currentValue)
        
    b = datetime.datetime.now()
    print(f"The program started at {b.strftime('%Y-%m-%d %H.%M')}")

    t_start = process_time()
    # parsing results
    result_handle = open(file)
    blast_records = NCBIXML.parse(result_handle)

    result_table =[]
    
    for alignment in blast_records:
        query=alignment.query # ID of query
        # checking of possible alignment
        if (len(alignment.alignments))!=0:
            for hit in alignment.alignments:
                for hsp in hit.hsps:
                    e = hsp.expect # e-value
                    # Filtering results for output
                    if e<evalue:
                        result_table.append([query,alignment.query_letters,hsp.align_length,hit.length,hit.title,e,hsp.score,hsp.align_length,hsp.identities])
                    else:
                       result_table.append([query,alignment.query_letters,hsp.align_length,hit.length,'No',e,hsp.score,hsp.align_length,hsp.identities]) 

    df = pd.DataFrame(columns=['ID_query','query_length','Aligments length','Hit length','Blast_hit', 'e_value','Score','align_len','identities'], data=result_table)

    df.to_csv(r'{}_table_{}.csv'.format(file,evalue), index=False, compression='gzip')
    t_stop = process_time()

    print(f"TIME: {(t_stop-t_start)/60} min") 


