from Bio.Blast import NCBIXML
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO



# path to BLAST results
out_file = r"C:\Users\fesenkoi2\IFESENKO\NCBI_BLASTX_TBLASTN\sORFs_filtered_lncRNAs_translation_nonredundant.fasta.out5tblastn_mosses"

physcomitrella_file = r"C:\Users\fesenkoi2\IFESENKO\orthologs\physcomitrella_ortho_mosses_tblastn.fasta"
physcomitium_file = r"C:\Users\fesenkoi2\IFESENKO\orthologs\ortho_mosses_tblastn.fasta"

records_physco = []
records_physcomitrium = []

cnt_all = set()
cnt_hom = set()

# parsing results
result_handle = open(out_file)
blast_records = NCBIXML.parse(result_handle) 
results={}
print('START...')
table = set()

    
for alignment in blast_records:
    query=alignment.query # ID of query
    cnt_all.add(query) # count of records
    
    # checking of possible alignment
    if (len(alignment.alignments))!=0:
        for hit in alignment.alignments:
            for hsp in hit.hsps:
                e = hsp.expect # e-value
                # Filtering results for output
                if e<0.00001:
                    cnt_hom.add(query)
                    h=hit.title
                    table.add((query,h))
                    query_aligm = Seq(hsp.query, IUPAC.unambiguous_dna)
                    hit_aligm  = Seq(hsp.sbjct, IUPAC.unambiguous_dna)
                    records_physco.append(SeqRecord(seq=query_aligm, id=query))
                    records_physcomitrium.append(SeqRecord(seq=hit_aligm, id=h))
                    cnt_hom.add(query)    
                      
                       
                        
SeqIO.write(records_physco,physcomitrella_file,"fasta")

SeqIO.write(records_physcomitrium,physcomitium_file,"fasta")

with open (r'C:\Users\fesenkoi2\IFESENKO\orthologs\sORFs_table_orthologs.txt', 'w') as out_file:
    for pair in table:
        out_file.write(f"{pair[0]} {pair[1]}\n")


print(f"The number of query: {len(cnt_all)}, passed the filter: {len(cnt_hom)}")