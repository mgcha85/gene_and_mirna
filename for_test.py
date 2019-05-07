from Database import Database
import os
import sqlite3
import pandas as pd

root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'

def split_by_tissues():
    fpath = os.path.join(root, 'database/Fantom/v5', 'hg19.fantom_cross_check_ucsc_ensembl_gencode_out.db')
    con = sqlite3.connect(fpath)
    tlist = Database.load_tableList(con)

    for tname in tlist:
        print(tname)
        df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)

        result = {}
        for idx in df.index:
            tissues = df.loc[idx, 'tissues'].split(';')
            start = df.loc[idx, 'fan_start']
            end = df.loc[idx, 'fan_end']
            resources = df.loc[idx, 'resources']

            for tis in tissues:
                if tis not in result:
                    result[tis] = [[start, end, resources]]
                else:
                    result[tis].append([start, end, resources])

        for tis, contents in result.items():
            con_out = sqlite3.connect(fpath.replace('.db', '_{}.db'.format(tis)))
            df_res = pd.DataFrame(data=contents, columns=['start', 'end', 'resources'])
            df_res.to_sql(tname, con_out, if_exists='replace', index=None)

split_by_tissues()