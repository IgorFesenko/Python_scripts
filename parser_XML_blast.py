from Bio.Blast import NCBIXML
import pandas as pd


# path to BLAST results
out_file = r"C:\Users\fesenkoi2\IFESENKO\SPADA\spada_transcripts_translation_blastp_physco"

# path to alignment file
out_aligm = r"C:\Users\fesenkoi2\IFESENKO\SPADA\spada_transcripts_translation_blastp_physco_aligment.txt"

# create table
df = pd.DataFrame({'ID_query':[], 'query_length':[],'Aligments length':[], 'Hit length':[], 'Blast_hit':[], 'e_value':[], 'Score':[]})

#counters
cnt_all = set()
cnt_hom = set()

# parsing results
result_handle = open(out_file)
blast_records = NCBIXML.parse(result_handle) 
results={}
print('START...')

with open(out_aligm, 'a') as out_file:
    
    for alignment in blast_records:
        query=alignment.query # ID of query
        cnt_all.add(query) # count of records
        l=alignment.query_letters # the length of query

        # checking of possible alignment
        if (len(alignment.alignments))!=0:
            
            for hit in alignment.alignments: 
                for hsp in hit.hsps:
                    p = hsp.identities # the number of identities
                    d = hit.length # hit length
                    c = hsp.align_length # alignments length
                    e = hsp.expect # e-value
                                       
                    # Filtering results for output
                    if e<0.00001 and p/l*100>=70: # and :p/d*100>=90 and c/l*100>=90
                        h=hit.title # Blast hit
                        score=hsp.score # Blast score
                        cnt_hom.add(query)
                      
                        # writing alignment to a file
                        out_file.write('****Alignment****'+'\n')
                        out_file.write('query: {}'.format(query)+'\n')
                        out_file.write('Blast_hit: {}'.format(h)+'\n')
                        out_file.write('e value: {}, Score: {}'.format(e, score)+'\n')
                        out_file.write('identities: {} length query: {} aligments length: {}'.format(p,l,c)+'\n')
                        out_file.write('Hit length: {}'.format(d)+'\n')
                        out_file.write(hsp.query[0:100] + '...'+'\n')
                        out_file.write(hsp.match[0:100] + '...'+'\n')
                        out_file.write(hsp.sbjct[0:100] + '...'+'\n')
                        out_file.write('\n')
                        
                        # add to output table
                        df = df.append({'ID_query':query, 'query_length':l,'Aligments length':c, 'Hit length':d,'Blast_hit':hit, 'e_value':e, 'Score':score}, ignore_index=True)
                        
                    else:
                        # add to output table
                        df = df.append({'ID_query':query, 'query_length':l,'Aligments length':c, 'Hit length':d,'Blast_hit':'No', 'e_value':'No', 'Score':'No'}, ignore_index=True)
        else:
            # add to output table
            df = df.append({'ID_query':query, 'query_length':l,'Aligments length':'No', 'Hit length':'No','Blast_hit':'No', 'e_value':'No', 'Score':'No'}, ignore_index=True)


# writing output table
df.to_excel(r"C:\Users\fesenkoi2\IFESENKO\SPADA\spada_transcripts_translation_blastp_physco_eval00001_align70.xlsx", index=False)

print(f"The number of query: {len(cnt_all)}, passed the filter: {len(cnt_hom)}")