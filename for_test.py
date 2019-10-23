import pandas as pd
import sqlite3

# fpath = '/home/mingyu/Downloads/miRNA_Tss_100_redefined_4.bed'
# df = pd.read_csv(fpath, names=['miRNA', 'loc'], sep='\t')
# df['loc'] = df['loc'].str[1:-1]
# df['loc'] = df['loc'].str.split(', ')
# for idx in df.index:
#     loc = df.loc[idx, 'loc']
#
#     ll = []
#     for l in loc:
#         num, chr, start, end = l.split('_')
#         ll.append('_'.join([chr, start, end]))
#     df.loc[idx, 'loc'] = ';'.join(ll)
#
# df.to_csv(fpath, sep='\t', index=None, header=False)
# exit(1)

# ref
fpath = '/home/mingyu/Bioinformatics/papers_tss.db'
con_mir = sqlite3.connect(fpath)
df_mir = pd.read_sql("SELECT * FROM 'PRO'", con_mir, index_col='name')

fpath = '/home/mingyu/Bioinformatics/database/mirbase/22/mirbase.db'
con_mir = sqlite3.connect(fpath)
df_mirbase = pd.read_sql("SELECT * FROM 'hsa_liftover_hg19' WHERE feature='miRNA_primary_transcript'", con_mir, index_col='Name')

fpath = '/home/mingyu/Bioinformatics/database/consistent_miRNA_330.bed'
df = pd.read_csv(fpath, sep='\t', names=['miRNA', 'loc'])
loc = df['loc'].str.split(';')
for idx in loc.index:
    mir = df.loc[idx, 'miRNA']
    row = loc[idx]
    for ele in row:
        chr, start, end = ele.split('_')
        df.loc[idx, 'chromosome'] = chr
        df.loc[idx, 'start'] = start
        df.loc[idx, 'end'] = end
        if mir in df_mir.index:
            df_sub = df_mir.loc[mir:mir, :]
            if df_sub.shape[0] > 0:
                df_sub = df_sub.iloc[:1, :]
            df.loc[idx, 'strand'] = df_sub.loc[mir, 'strand']
        else:
            df.loc[idx, 'strand'] = df_mirbase.loc[mir, 'strand']

con_out = sqlite3.connect(fpath.replace('.bed', '.db'))
for chr, df_chr in df.groupby('chromosome'):
    for str, df_str in df_chr.groupby('strand'):
        df_str.drop('loc', axis=1).to_sql('_'.join([chr, str]), con_out, index=None)
