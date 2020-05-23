import requests
from bs4 import BeautifulSoup
import re


main_url = r'https://data.cyverse.org/dav-anon/iplant/projects/commons_repo/curated/oneKP_capstone_2019/transcript_assemblies'

r = requests.get(main_url)
soup1 = BeautifulSoup(r.text, 'html.parser')
directories = [main_url + '/' + node.get('href') for node in soup1.find_all('a')]

# reading the list of species
with open (r"C:\Users\fesenkoi2\IFESENKO\1000_transcriptomes\Conifers.txt") as input_lst:
    names = [i.strip().split('-')[0] for i in input_lst.readlines()]

print(f'The number of species:{len(names)}')
print()


for d in directories: #reading directories
    n = re.findall(r'[A-Z]{4}',d)
    if len(n)>0 and n[0] in names:
        print(d)   # reading list of names
        req = requests.get(d)
        print(req)
        try:
            soup = BeautifulSoup(req.text, 'html.parser')
            files = [d + '/' + node.get('href') for node in soup.find_all('a')]
            for i in files:
                if re.findall('Trans-assembly', i):
                    myfile = requests.get(i)
                    name = d.split('/')[-2]
                    print(name)
                    print()
                    open (r'C:\Users\fesenkoi2\IFESENKO\1000_transcriptomes\Conifers\{}.fa.bz2'.format(name), 'wb').write(myfile.content)
        except:
            print('ERROR:',d)
            




