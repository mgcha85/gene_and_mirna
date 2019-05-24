from Database import Database
import os
import sqlite3
import pandas as pd
import numpy as np


fpath = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/software/stringtie-1.3.3b/genes.gtf'
gtf_columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

df = pd.read_csv(fpath, names=gtf_columns, comment='#', sep='\t')
df['chromosome'] = 'chr' + df['chromosome'].astype(str)
attribute = df['attribute'].str.split('; ')

con = sqlite3.connect('/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/software/stringtie-1.3.3b/genes.db')

pkm = []
for attr in attribute:
    dict = {}
    for a in attr[:-1]:
        key, value = a.split(' ')
        dict[key] = value.replace('"', '')

    if 'gene_name' not in dict:
        pkm.append(np.nan)
        continue
    pkm.append(dict['gene_name'])

df_res = pd.DataFrame(data=pkm, columns=['gene_name'])
df = pd.concat([df, df_res], axis=1)
df.drop(['attribute', 'frame'], axis=1).to_sql('gencode.v30lift37.basic.annotation', con, index=None, if_exists='replace')
