import pandas as pd
import sqlite3
import socket
import os


class high_consistent_genes:
    def __init__(self, root):
        self.root = root

    def to_server(self, root, tag):
        from Server import Server
        import sys

        which = 'newton'
        server = Server(root, which=which)
        server.connect()

        local_path = sys.argv[0]
        dirname, fname__ = os.path.split(local_path)
        fname = fname__.replace('.py', '_{}.py'.format(tag))

        curdir = os.getcwd().split(os.sep)[-1]
        server_root = os.path.join(server.server, 'source', curdir).replace(os.sep, '/')
        server_path = local_path.replace(dirname, server_root)
        server_path = server_path.replace(fname__, fname)

        server.job_script(fname, src_root=server_root, time='08:00:00', pversion=3)
        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', server_root + '/' + 'dl-submit.slurm')

        stdin, stdout, stderr = server.ssh.exec_command("cd {root};sbatch {root}/dl-submit.slurm".format(root=server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def run(self):
        fpath = os.path.join(self.root, 'database/Fantom/v5/tissues', 'correlation_fan_rna_100_processed.xlsx')
        df = pd.read_excel(fpath, index_col=0)
        df = df[df['corr'] >= 0.75]

        df_res = []
        for gene, df_grp in df.groupby('gene_name'):
            idx = df_grp['corr'].idxmax()
            df_res.append(df_grp.loc[idx])
        df_res = pd.concat(df_res, axis=1).T
        df_res.to_excel(fpath.replace('.xlsx', '_by_gene.xlsx'))


if __name__ == '__main__':
    hostname = socket.gethostname()
    if hostname == 'mingyu-Precision-Tower-7810':
        root = '/home/mingyu/Bioinformatics'
    elif hostname == 'DESKTOP-DLOOJR6' or hostname == 'DESKTOP-1NLOLK4':
        root = 'D:/Bioinformatics'
    elif hostname == 'mingyu-Inspiron-7559':
        root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
    else:
        root = '/lustre/fs0/home/mcha/Bioinformatics'

    hcg = high_consistent_genes(root)
    if hostname == 'DESKTOP-DLOOJR6' or hostname == '-1NLOLK4':
        hcg.to_server(root, "")
    else:
        hcg.run()
