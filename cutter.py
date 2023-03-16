# Cut fasta db on fragments

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# read file
file_input = r""

SEQ_LIMIT = 990

names = set()
def fasta_record(rec,cnt):
    name = '{}_fragment_{}.fasta'.format(file_input.split('.')[0],cnt)
    names.add(name)
    print(name)
    SeqIO.write(rec,name,"fasta")

seq_count=0
records=[]
cnt=0

for record in SeqIO.parse(file_input, 'fasta'):
    cnt+=1
    seq_count+=1
    if seq_count<SEQ_LIMIT:
        records.append(SeqRecord(seq=record.seq, id=record.id, description=record.description))
            
    else:
        fasta_record(records,cnt)
        records=[]
        seq_count=0

fasta_record(records,cnt)
print('The number of records: ',cnt)
