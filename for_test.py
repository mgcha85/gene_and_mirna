import pandas as pd
import sqlite3

df1 = pd.read_excel('/home/mingyu/Bioinformatics/database/gencode/high_correlated_fan_rna_100.xlsx')
df1['transcript_id'] = df1['transcript_id'].str.split('_').str[0]
df1 = df1.set_index('transcript_id', drop=True)
df2 = pd.read_csv('/home/mingyu/Downloads/correlated_transcripts.txt', sep='\t', index_col=0)

inter = set.intersection(set(df1.index), set(df2.index))
print(inter)