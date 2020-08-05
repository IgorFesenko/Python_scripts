#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Analysing alignmets after TBLASTN (web form)
'''
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
#import swifter
from time import process_time
import re



def benchmark(func):
    """
    Декоратор выводящий время работы функции
    """
    
    def wrapper(*args, **kwargs):
        t_start = process_time()
        res = func(*args, **kwargs)
        t_stop = process_time()
        print(f"{func.__name__}Time: {t_stop-t_start} sec\n")
        return res
    return wrapper

@benchmark
def read_seq(file, e):
    """
    Парсим файл вывода после 1Kb
    """
    blast_qresult = SearchIO.read(file, "blast-text")
    results = []
    for hit in blast_qresult:
        for hsp in hit:
            if hsp.evalue<e:
                seq_id = hit.id.split(':')[1]
                if int(hsp.query_frame)<0:
                    desc = 'reverse'
                else:
                    desc = 'true'
                results.append(SeqRecord(seq=hsp.hit.seq, id=f"gnl|onekp|{seq_id}", description=f"{desc}_{hit.description}"))
                #sequence = Seq(hsp.hit.seq, IUPAC.unambiguous_dna)           
    return results

@benchmark
def seq_table(file):
    """
    reading fasta file with results
    """
    
    seq_data = []
    for record in SeqIO.parse(file,'fasta'):
        descr = record.description.split()[1]
        descr1 = descr.split('_',maxsplit=1)[0]
        descr2 = descr.split('_',maxsplit=1)[1]
        seq_data.append([record.id, str(record.seq), descr1, descr2])
    return seq_data


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
                    
def longest_taq(row):
    # creating a tag for seqrching
    lst = row['Sequence'].split(r'-')
    return sorted(lst, key=len, reverse=True)[0]


def transcript_extraction(row):
    # adding transcript for fasta
    seq_seq = dict_trascript[row['ID']].seq
    #seq_seq = Seq(seq_transcipt[start:stop], IUPAC.unambiguous_dna)
    if row['frame'] == 'reverse':
            return str(seq_seq.reverse_complement())
    else:
        return str(seq_seq)
    
def transcr_refine(row):
    if 'N' in str(row['transcript']):
        cng_seq = str(re.sub('N','A', str(row['transcript'])))
        return cng_seq
    return str(row['transcript'])


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
    e = float()
        
       
    for currentArgument, currentValue in arguments:
        if currentArgument in ("-d", "--dir"):
            directory = currentValue                     
        elif currentArgument in ("-e", "--eval"):
            e = float(currentValue)
        


    E_VALUE = e

    print(f"Directory: {directory}")

    # tblastn file
    file = f"{directory}/out.0.txt"

    table_out = f"{directory}/script_results/table_seq.xlsx"

    # all transcript sequences
    transcripts = f"{directory}/all_seq.fasta"
    dict_trascript = SeqIO.index(transcripts, "fasta")

    #creating directory
    
    os.system(f"mkdir {directory}/script_results")


    print("reading search results")
    res = read_seq(file,E_VALUE)


    # wrting fasta results
    seq_name = f"{directory}/fasta_aligments.fasta"
    SeqIO.write(res,seq_name,"fasta")

    print("reading fasta")
    seq = seq_table(seq_name)

    # creating table
    df = pd.DataFrame(columns=['ID','Sequence','frame', 'species'], data=seq)
    df['taq'] = df.apply(longest_taq, axis=1)


    # adding transcript
    df['transcript'] = df.apply(transcript_extraction, axis=1)

    #refine transcript
    df['transcript'] = df.apply(transcr_refine, axis=1)


    #adding full peptide sequence
    df['Full_PEP'] = df.apply(transcript_translation, axis=1)

    print(df[['ID','taq','frame','species','Full_PEP']])

    df.to_excel(table_out, index=False)

    # reading table
    filt_table = df[np.logical_not(df['Full_PEP']==None)]
    out_fasta = f"{directory}/script_results/full_seq.fasta"

    results=[]

    #reading peptide sequence and name
    with open (f"{directory}/info.txt") as inp:
        lst_pep = [i.strip() for i in inp.readlines()]
    seq_seq = Seq(lst_pep[1], IUPAC.unambiguous_dna)
    results.append(SeqRecord(seq=seq_seq, id=lst_pep[0]))

    #creating fasta
    for index,row in filt_table.iterrows():
        seq_id = f"{row['ID']}_{row['species']}"
        if row['Full_PEP'] != 'None':
            seq_seq = Seq(row['Full_PEP'], IUPAC.unambiguous_dna)
            results.append(SeqRecord(seq=seq_seq, id=seq_id))
        else:
            continue
            

        

    SeqIO.write(results,out_fasta,"fasta")



        
    
    