""""
Удаляем * в файлах баз данных
"""
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from time import process_time

t_start = process_time()

#file with initial sequences 
file_input = r'/Users/igorfesenko/rab_proteins.fasta'

out_file = r'/Users/igorfesenko/rab_proteins_asterix.fasta'

records = []

for record in SeqIO.parse(file_input, 'fasta'):
    sequence = Seq(str(record.seq).replace('*',''))
    records.append(SeqRecord(seq=sequence, id=record.id, description=record.description))

SeqIO.write(records, out_file,"fasta")


t_stop = process_time()
print(f"Time: {t_stop-t_start} sec")