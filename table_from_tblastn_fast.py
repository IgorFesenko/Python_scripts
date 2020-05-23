'''
Creating a table with sequences and ID for orthologs
'''

import pandas as pd
from Bio.Blast import NCBIXML
from time import process_time
import datetime

b = datetime.datetime.now()
print(f"The program started at {b.strftime('%Y-%m-%d %H.%M')}")

# path to BLAST results
out_file = r"C:\Users\fesenkoi2\IFESENKO\data_for_orthologs\filt_proteincoding_proteins.fasta.out5tblastn_mosses"


result_handle = open(out_file)
blast_records = NCBIXML.parse(result_handle) 

print('START...')
result_table=[]
t_start = process_time()
    
for alignment in blast_records:
    query=alignment.query # ID of query
    # checking of possible alignment
    if (len(alignment.alignments))!=0:
        for hit in alignment.alignments:
            for hsp in hit.hsps:
                e = hsp.expect # e-value
                # Filtering results for output
                if e<0.001:
                    result_table.append([query,hit.title,hsp.query,hsp.sbjct,'-','-',e,hsp.query_start, hsp.query_end,hsp.sbjct_start,hsp.sbjct_end, hsp.frame, hsp.align_length,hsp.identities]) 

# writing first output table
print('WRITING first output...')
df = pd.DataFrame(columns=['query','hit','Seq_query','Seq_hit', 'Seq_query_nucl','Seq_hit_nucl','evalue','query_start','query_stop','hit_start','hit_stop', 'frame','align_len','identities'], data=result_table)
print(df.head()) 

df.to_csv(r'C:\Users\fesenkoi2\IFESENKO\data_for_orthologs\filt_proteincoding_proteins_tblastn_mosses_table_10e3.csv', index=False, compression='gzip')
t_stop = process_time()

print(f"TIME: {(t_stop-t_start)/60} min")