#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Анализируем данные TBLASTN полученные в web форме
"""

import os
import sys
import getopt
from Bio import SearchIO
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
from Bio.Blast import NCBIXML
from time import process_time
import re

def longest_taq(row):
    # creating a tag for searching
    lst = row['Sequence'].split(r'-')
    if lst:
        return sorted(lst, key=len, reverse=True)[0]
    else:
       return None

def transcript_extraction(row):
    # adding transcript for fasta
    seq_seq = dict_trascript[row['ID']].seq
    #seq_seq = Seq(seq_transcipt[start:stop], IUPAC.unambiguous_dna)
    if float(row['frame']) < 0:
            return str(seq_seq.reverse_complement())
    else:
        return str(seq_seq)

def transcr_refine(row):
    if 'N' in str(row['transcript']):
        cng_seq = str(re.sub('N','A', str(row['transcript'])))
        return cng_seq
    return str(row['transcript'])

def transcript_translation(row):
    """
    Transcript translation in 3-frame
    """
    sequence = Seq(row['transcript'], IUPAC.unambiguous_dna)
    
    translation1 = sequence.translate()
    #print(translation1)
    translation2 = sequence[1:].translate()
    translation3 = sequence[2:].translate()
    lst_translation = [translation1,translation2,translation3]
    peptide = pep_in_translation(lst_translation,row['taq'])
    #print(type(peptide))
    return str(peptide)


def pep_in_translation(pep_lst, taq):
    
    """
    Finding peptide in translation list
    """
    
    for i in range(3):
        for pep in pep_lst[i].split('*'):
            if 'M' in pep:
                position = pep.find('M')
                new_seq = pep[position:]
                #print(new_seq)
                if taq in new_seq:
                    return new_seq
    return None 

if __name__ == "__main__":
    
    unixOptions = "d:e:"  
    gnuOptions = ["dir=",'eval=']
   
    fullCmdArguments = sys.argv
    argumentList = fullCmdArguments[1:]

   
    try:  
        arguments, values = getopt.getopt(argumentList, unixOptions, gnuOptions)
    except getopt.error as err:  
        print (str(err))
        sys.exit(2) 

    
    directory = ''
    evalue = float()
        
       
    for currentArgument, currentValue in arguments:
        if currentArgument in ("-d", "--dir"):
            directory = currentValue                     
        elif currentArgument in ("-e", "--eval"):
            evalue = float(currentValue)

    print('Reading transcripts...')
    # all transcript sequences
    transcripts = f"{directory}/all_seq.fasta"
    dict_trascript = SeqIO.index(transcripts, "fasta")

    
    print('Start parsing XML...')

    t_start = process_time()

    file = f"{directory}/out.5.xml"

    result_handle = open(file)
    blast_records = NCBIXML.parse(result_handle)

    
    
    for alignment in blast_records:
        qry = alignment.query # ID of query
        gry_name = qry.split(' ')[0]
        #replace characters
        gry_name = gry_name.replace(".", '_')
        gry_name = gry_name.replace("|", '_')
        gry_name = gry_name.replace(":", '_')
        
        # checking of possible alignment
        if len(alignment.alignments) != 0:
            result_table =[]
            
            for hit in alignment.alignments:
                for hsp in hit.hsps:
                    e = hsp.expect # e-value
                    # Filtering results for output
                    if e < evalue:
                        os.system(f"mkdir {directory}/{gry_name}")    
                        result_table.append([qry.split(' ')[1],hit.title.split(' ')[0], hsp.sbjct, hsp.frame[1],hit.title.split(' ')[1]])                       
                        # creating table
                        print('Creating table...')
                        df = pd.DataFrame(columns=['query','ID','Sequence','frame', 'species'], data=result_table)
                        df['query'] = df['query'].map(lambda x: x.rsplit('_', maxsplit=1)[0])
                        #print(df.head())
                        df['taq'] = df.apply(longest_taq, axis=1)

                        # adding transcript
                        df['transcript'] = df.apply(transcript_extraction, axis=1)

                        #refine transcript
                        df['transcript'] = df.apply(transcr_refine, axis=1)

                        #adding full peptide sequence
                        df['Full_PEP'] = df.apply(transcript_translation, axis=1)

                        print(qry)

                        #print(df.head(5).to_string())
                        table_out = f"{directory}/{gry_name}/table_seq.xlsx"
                        df.to_excel(table_out, index=False)
                        print()

                                        
                                        
                        t_stop = process_time()

                        print(f"TIME: {(t_stop-t_start)/60} min") 






