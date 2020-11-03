# Extracting first exon-utr sequences from gff and writing in bed file

import pandas as pd

bed_file =[]
lst_utrs=[]

def utr_extract(df):
    for i in range(len(df)-1):
        name = df.loc[i,'feature']
        next_name = df.loc[i+1,'feature']
        if name=='five_prime_UTR':
            bed_file.append((df.loc[i,'seqname'],df.loc[i,'start'],df.loc[i,'end'],df.loc[i,'attribute'],'1000',df.loc[i,'strand']))
            if next_name=='CDS':
                bed_file.append((df.loc[i+1,'seqname'],df.loc[i+1,'start'],df.loc[i+1,'end'],df.loc[i+1,'attribute'],'1000',df.loc[i+1,'strand']))
        


input_file = r'/Users/igorfesenko/Google Диск/Phytozome_Ppatens_29012020/PhytozomeV12/Ppatens/annotation/Ppatens_318_v3.3.gene.gff3'

df = pd.read_csv(input_file, sep='\t', header=2, names=['seqname','source','feature', 'start','end','score','strand','frame','attribute'])
print(df.head())

utr_extract(df)

with open('/Users/igorfesenko/utr_exon.bed','w') as out:
    for rec in bed_file:
        out.write(f"{rec[0]}\t{rec[1]}\t{rec[2]}\t{rec[3]}\t{rec[4]}\t{rec[5]}\n")