# Extracting intron sequences from gff and read in bed file

import pandas as pd

bed_file =[]
lst_introns=[]


def intron_extraction2(df):
    for i in range(len(df)):
        name = df.loc[i,'feature']
        if name=='five_prime_UTR' or name=='mRNA':
            lst_introns=[]
        elif name=='CDS':
            lst_introns.append((df.loc[i,'seqname'],df.loc[i,'start'],df.loc[i,'end']))
        
        elif name=='three_prime_UTR':
            if df.loc[i,'strand']=='-':
                #print('reverse')
                bed_file.append(lst_introns[::-1])
            else:
                bed_file.append(lst_introns)
            
    



input_file = r"C:\Users\fesenkoi2\IFESENKO\sORFinder\Sfallax\Sfallax_522_v1.1.gene.gff3"

df = pd.read_csv(input_file, sep='\t', header=2, names=['seqname','source','feature', 'start','end','score','strand','frame','attribute'])

df.head()
intron_extraction2(df)


for gene in bed_file:
    if len(gene) > 1:
        for i in range(len(gene)-1):
            try:
                chr = gene[i][0]
                intr_start = gene[i][2]+1
                intr_end = gene[i+1][1]-1
                with open(r'C:\Users\fesenkoi2\IFESENKO\sORFinder\Sfallax\sfallax_intron2.bed','a') as out:
                    out.write(f"{chr}\t{intr_start}\t{intr_end}\n")
            except:
                print(gene)
                continue


                    