#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Parsing HMMSEARCH output
'''
from Bio import SearchIO
import pandas as pd
import sys
import getopt

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
        
    
    # parsing results
    results = []
    for qresult in SearchIO.parse(file, "hmmer3-text"):
        for hit in qresult:
            for hsp in hit:
                if hit.evalue<evalue:
                    results.append([qresult.id,str(hsp.query.seq),hit.id,hit.description,str(hsp.hit.seq),hsp.aln_span,hit.evalue, hsp.evalue, hsp.bitscore])
                
    # table creating
    df = pd.DataFrame(columns=['Query','Query_seq','Hit_id', 'Hit_description','Hit_seq','aln_length','hit_evalue','hsp_evalue','Bit_score'], data=results)

    #output
    df.to_csv(f'./{file}.xlsx', index=False)