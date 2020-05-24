'''
Parsing the results of TMHMM2 output
'''
import pandas as pd
import re

file = r"C:\Users\fesenkoi2\IFESENKO\TMHMM\NcbiFinder_output"

table = pd.DataFrame({'sORF':[],'Number_of_predicted_THMs':[],'Exp_number_AA':[],'Exp_number_first60AA':[],'Prob_of_N':[],'Possible_N-term_signal':[],'Result':[]})

with open(file) as inp_file:
    content = [i.strip() for i in inp_file.readlines()]

print(content[:10])
overall=[]

for line in content:
    if re.findall('Length',line):
        pred=""
        name = line.split(' ')[1]
        length = line.split(' ')[-1]
        
    elif re.findall('POSSIBLE', line):
        pred = re.split(pattern=r'ORF_\d+ ', string=line)[-1]

    elif line.startswith('O'):
        #print(line)
        domain = line.split('\t')[-2]
        span = line.split('\t')[-1]
        overall.append([name,length, domain,span,pred])

print(overall[:20])


table = pd.DataFrame(columns= ['sORF', 'Length','position','span','prediction'], data=overall)
print(table.head())
table.to_csv(r'C:\Users\fesenkoi2\IFESENKO\TMHMM\NcbiFinder_output.csv', index=False)



'''
for i in overall[:20]:
    if len(i)>7:
        print(i)
for line in content:
    if 'Length' in line:
        lst = []
        lst.append(line.split(' ')[1])
        lst.append(line.split(' ')[-1])

    elif 'N-term' in line:
           lst.append(line.split(' ')[2]) 
    elif 'TMHMM2' in line:
        lst.append(line.split('\t')[2])
        lst.append(line.split('\t')[3])

    else:
        lst.append(line.split(' ')[-1])
'''