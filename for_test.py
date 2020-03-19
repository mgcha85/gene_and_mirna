import pandas as pd
import sqlite3
from Database import Database
import os
from collections import OrderedDict


dirname = '/home/mingyu/Bioinformatics/database/RNA-seq/fastq'
flist = os.listdir(dirname)
fids = [os.path.splitext(x)[0] for x in flist]

def get_pair(contents):
    pair_list = OrderedDict()
    for url in contents:
        key = url.split('/')[-2]
        if key not in pair_list:
            pair_list[key] = [url]
        else:
            pair_list[key].append(url)
    return pair_list

with open('RNA-seq_list (copy).csv', 'rt') as f:
    rna_seq_all = f.read().split('\n')

with open('RNA-seq_list.csv', 'rt') as f:
    rna_seq_2 = f.read().split('\n')

fids = set(fids) #set.union(set(get_pair(rna_seq_2).keys()), set(fids))

result = []
for key, value in get_pair(rna_seq_all).items():
    if key not in fids:
        result += value

with open('RNA-seq_list2.csv', 'wt') as f:
    f.write('\n'.join(sorted(result)))
exit(1)


# fpath = '/home/mingyu/Bioinformatics/database/gencode/high_correlated_fan_rna_100__.db'
# con = sqlite3.connect(fpath)
# dfs = []
# for tname in Database.load_tableList(con):
#     dfs.append(pd.read_sql("SELECT * FROM '{}'".format(tname), con))
# pd.concat(dfs).to_excel(fpath.replace('.db', '.xlsx'))
# exit(1)

# fpath = '/home/mingyu/Bioinformatics/database/Fantom/v5/tissues/correlation_fan_rna_100_raw.db'
# con = sqlite3.connect(fpath)
#
# dfs = []
# for tname in Database.load_tableList(con):
#     dfs.append(pd.read_sql("SELECT * FROM '{}'".format(tname), con))
# df = pd.concat(dfs)
# print(df.shape[0])
# exit(1)

fpath = "/home/mingyu/Bioinformatics/database/gencode/high_correlated_fan_rna_100.xlsx"
df = pd.read_excel(fpath)
gene_names1 = sorted(list(set(df['gene_name'])))

mean = df['corr'].mean()
median = df['corr'].median()
std = df['corr'].std()
min = df['corr'].min()
max = df['corr'].max()
q7 = df['corr'].quantile(0.7)
q8 = df['corr'].quantile(0.8)
q9 = df['corr'].quantile(0.9)

print('mean: {:0.2f}, median: {:0.2f}, std: {:0.2f}, min: {:0.2f}, max: {:0.2f}, q7: {:0.2f}, q8: {:0.2f}, q9:{:0.2f}'
      ''.format(mean, median, std, min, max, q7, q8, q9))

with open('/home/mingyu/Bioinformatics/database/gencode/high_correlated_genes.txt', 'rt') as f:
    gene_names2 = sorted(list(set(f.read().split('\n'))))

inter = set.intersection(set(gene_names1), set(gene_names2))
print('1: {}, 2: {}, inter: {}'.format(len(gene_names1), len(gene_names2), len(inter)))
