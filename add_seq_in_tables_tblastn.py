'''
Adding nucleotide sequences in table
'''
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
from time import process_time
import datetime


script = 'table_from_tblastn_fast.py'

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
        return str(seq_seq.reverse_complement())
    else:
        return seq_seq

print('READING TABLE...')

table_file=r"C:\Users\fesenkoi2\IFESENKO\data_for_orthologs\vill_mipepid_tblastnAUG_nostop_name.csv"

df = pd.read_csv(table_file, compression='gzip')


print('READ DATABASE...')
t_start2 = process_time()

db1 = r"C:\Users\fesenkoi2\IFESENKO\data_for_orthologs\Villersexel.fa"
db2 = r"C:\Users\fesenkoi2\IFESENKO\lncRNAs_sORFs\lncRNAs_sorfs_mipepid\sORFs_mipepid_nucleotide_db_newname.fasta"


dict_db1 = SeqIO.index(db1, "fasta") # mosses transcripts
dict_db2 = SeqIO.index(db2, "fasta") # mipepid


print("ADD NUCLEOTIDES...")
print('adding queries..')
df['query_nucleotide'] = df.apply(query_sequence,axis=1)

print('adding hits...')
df['hit_nucleotide'] = df.apply(hit_sequence,axis=1)


df.to_csv(r"C:\Users\fesenkoi2\IFESENKO\data_for_orthologs\table_sorfs_vill_tblastnAUG_nostop.csv_nucl.csv", compression='gzip', index=False)


'''
records_physco = []
records_mosses = []

cnt=0
for rec in table:
    
        #print(rec[0])
    seq_transcript = str(dict_db1[rec[1]].seq)
    seq_sorf = str(dict_db2[rec[0]].seq)
        #print(seq_seq)
    seq_seq = seq_transcript[int(rec[2])-1:int(rec[3])]

    query_aligm = Seq(seq_sorf, IUPAC.unambiguous_dna)
    hit_aligm  = Seq(seq_seq, IUPAC.unambiguous_dna)
    records_physco.append(SeqRecord(seq=query_aligm, id=rec[0]))
    records_mosses.append(SeqRecord(seq=hit_aligm, id=rec[1]))

    seq_transcript = str(dict_db1[rec[1]].seq)
    seq_sorf = str(dict_db2[rec[0]].seq)
    seq_seq = seq_transcript[int(rec[2])-1:int(rec[3])]
    #sequences
    query_aligm = Seq(seq_sorf, IUPAC.unambiguous_dna)
    hit_aligm  = Seq(seq_seq, IUPAC.unambiguous_dna)

    if rec[4]<0:
        records_mosses.append(SeqRecord(seq=hit_aligm.reverse_complement(), id=rec[1]))
    else:
        records_mosses.append(SeqRecord(seq=hit_aligm, id=rec[1]))
    records_physco.append(SeqRecord(seq=query_aligm, id=rec[0]))
        
'''

t_stop2 = process_time()
print(f"TIME: {(t_stop2-t_start2)/60} min")
