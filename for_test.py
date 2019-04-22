# import pandas as pd
# import numpy as np
# from scipy.stats import spearmanr
#
# n_rows = 2500
# cols = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
#
# df = pd.DataFrame(np.random.random(size=(n_rows, len(cols))), columns=cols)
# v = np.random.random(size=len(cols))
#
# # original implementation
# corr, _ = zip(*df.apply(lambda x: spearmanr(x,v), axis=1))
# corr = pd.Series(corr)
#
# # modified implementation
# df1 = df.rank(axis=1)
# v1 = pd.Series(v, index=df.columns).rank()
# corr1 = df1.corrwith(v1, axis=1)
#
# print(corr1)

import sqlite3
import pandas as pd
import os

# dirname = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/database/ensembl/TSS'
# columns = ['Chromosome/scaffold name', 'Transcript start (bp)', 'Transcript end (bp)']
#
# con = sqlite3.connect(os.path.join(dirname, 'mart_export.db'))
# df = pd.read_sql_query("SELECT * FROM 'Ensembl'", con)
#
# df_hg19 = pd.read_csv(os.path.join(dirname, 'Ensembl_hg19.bed'), sep='\t', names=columns)
#
# df[columns] = df_hg19[columns]
# df_str = df.groupby('Strand')
#
# dfs = []
# for str, df_sub in df_str:
#     if str > 0:
#         df_sub['Transcription start site (TSS)'] = df_sub['Transcript start (bp)']
#     else:
#         df_sub['Transcription start site (TSS)'] = df_sub['Transcript end (bp)']
#     dfs.append(df_sub)
#
# df = pd.concat(dfs)
# df['Chromosome/scaffold name'] = 'chr' + df['Chromosome/scaffold name']
#
# out_con = sqlite3.connect(os.path.join(dirname, 'mart_export2.db'))
# df.to_sql('Ensembl', out_con, index=None)

df = pd.read_csv('/home/mingyu/Downloads/FPKM_allsamples.txt', sep='\t')

