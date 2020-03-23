import sqlite3
import pandas as pd
import os
#
# root = 'D:/Bioinformatics/database/Fantom/v5/tissues'
# fpath = os.path.join(root, 'correlated_transcripts.txt')
# df_aml = pd.read_csv(fpath, sep='\t', names=['gid', 'tid', 'corr'], index_col=1)
# df_min = pd.read_excel(os.path.join(root, 'correlation_fan_rna_100_processed.xlsx'), index_col=0)
# df_min = df_min[df_min['corr'] > 0.75]
#
# A = set(df_min.index)
# B = set(df_aml.index)
# inter = set.intersection(A, B)
# remain = list(A - B)
# print(remain)
# remain = list(set(df_aml.index) - set(df_min.index))
# with open('remain.txt', 'wt') as f:
#     f.write('\n'.join(remain))

# exit(1)

# cols = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
# df = pd.read_csv("D:/Bioinformatics/database/gencode/gencode.v33lift37.annotation.gtf", sep='\t', comment='#', names=cols, chunksize=1<<4)
# for df_sample in df:
#     df_sample.to_csv("D:/Bioinformatics/database/gencode/gencode.v33lift37.annotation_sample.gtf", sep='\t')
#     exit(1)

with open('remain.txt', 'rt') as f:
    remains = f.read().split('\n')

fpath = "D:/Bioinformatics/database/gencode/gencode.v32lift37.annotation_attr_merged.db"
con = sqlite3.connect(fpath)
df = pd.read_sql("SELECT transcript_id, gene_type, transcript_type, feature FROM 'gencode.v32lift37' WHERE feature='transcript'", con, index_col='transcript_id')
# inter = set.intersection(set(remains), set(df.index))

A = set(remains)
B = set(df.index)
not_in_gencode = A - B
print(not_in_gencode)
inter = set.intersection(A, B)

df = df.loc[inter]
df_corr = pd.read_excel("D:/Bioinformatics/database/Fantom/v5/tissues/correlation_fan_rna_100_processed.xlsx", index_col=0)
for idx in df.index:
    if idx in df_corr.index:
        df.loc[idx, 'corr'] = df_corr.loc[idx, 'corr']

df.to_excel('temp.xlsx')
