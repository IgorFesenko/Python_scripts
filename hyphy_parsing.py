#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Parsing of json files from HyPhy (busted)

import json
import os
import pandas as pd
 

file_json = '/Users/igorfesenko/pal2nal_tables/tmp_fasta/mafft_align/NUCL_sORF_23229.phy.BUSTED.json'


if __name__ == "__main__":

    # results list
    results = []

    # reading list of json files
    names = [n for n in os.listdir() if n.endswith('json')]
    print(f"ANALYSIS : {len(names)} files")



    # Parsing files
    for n in names:
        with open (n) as f:
            parsed_json = json.load(f)
        res_data = parsed_json['test results'] # results {'LRT': 0, 'p-value': 0.5}
        inp_file = parsed_json['input']
        smorf = inp_file['file name'].split('NUCL_')[1].split('.')[0]
        results.append([smorf, res_data['LRT'], res_data['p-value']])

    df = pd.DataFrame(columns=['sORF','LRT', 'p_val'], data=results)

    df.to_csv(r'table_hyphy_from_json.csv', index=False)



