import sqlite3
import pandas as pd
import os
import socket
import zipfile


class Gencode:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def to_server(self):
        from Server import Server
        import sys

        which = 'newton'
        server = Server(self.root, which=which)
        server.connect()

        local_path = sys.argv[0]
        dirname, fname = os.path.split(local_path)
        curdir = os.getcwd().split('/')[-1]
        server_root = os.path.join(server.server, 'source', curdir)
        server_path = local_path.replace(dirname, server_root)

        server.job_script(fname, src_root=server_root, time='04:00:00')
        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def get_attr(self, attr):
        dfs = []
        for att in attr:
            at = [a.split(' ') for a in att]
            df = pd.DataFrame(at).set_index(0)
            df.iloc[:, -1] = df.iloc[:, -1].str.replace(';', '')
            df = df.loc[~df.index.duplicated(keep='first')].T
            dfs.append(df)
        return pd.concat(dfs).set_index(attr.index)

    def gtf_to_db(self, fpath, con):
        columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

        # chunksize = 1 << 20
        df_chunk = pd.read_csv(fpath, sep='\t', compression='gzip', names=columns, comment='#')
        # for df_chunk in pd.read_csv(fpath, sep='\t', compression='gzip', names=columns, chunksize=chunksize, comment='#'):
        attr = df_chunk['attribute'].str.replace('"', '').str.split('; ')
        df_attr = self.get_attr(attr)
        df = pd.concat([df_chunk.drop('attribute', axis=1), df_attr], axis=1)

        for chr, df_chr in df.groupby('chromosome'):
            for str, df_str in df_chr.groupby('strand'):
                tname = '_'.join([chr, str])
                try:
                    df_str.to_sql(tname, con, if_exists='append', index=None)
                except:
                    df_org = pd.read_sql("SELECT * FROM '{}'".format(tname), con)
                    df_str = pd.concat([df_str, df_org])
                    df_str.to_sql(tname, con, if_exists='replace', index=None)


if __name__ == '__main__':
    gc = Gencode()
    if gc.hostname == 'mingyu-Precision-Tower-7810':
        gc.to_server()

    else:
        fpath = os.path.join(gc.root, 'database/gencode/gencode.v32lift37.annotation.gtf.gz')
        con = sqlite3.connect(fpath.replace('.gtf.gz', '_2.db'))
        gc.gtf_to_db(fpath, con)
