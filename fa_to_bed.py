import os
import subprocess
import socket
import pandas as pd
import sqlite3
import numpy as np
import shutil
import sys
# from Server import Server


class fa2bed:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/home/mingyu/Bioinfomatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
            # self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.bowtie_root = os.path.join(self.root, 'software/hisat2-2.1.0')
        self.gtf_columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        self.bed_columns = ['chromosome', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'sequence', 'alignment', 'dummy1', 'dummy2', 'dummy3', 'dummy4', 'dummy5', 'dummy6', 'dummy7', 'dummy8', 'dummy9', 'dummy10', 'dummy11', 'dummy12']

    # def to_server(self):
    #     server = Server()
    #     server.connect()
    #
    #     local_path = sys.argv[0]
    #     dirname, fname = os.path.split(local_path)
    #
    #     server.job_script(fname, time='00:30:00')
    #
    #     server_root = os.path.join(server.server, 'source/gene_and_mirna')
    #     server_path = local_path.replace(dirname, server_root)
    #
    #     server.upload(local_path, server_path)
    #     server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))
    #
    #     stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
    #     job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
    #     print('job ID: {}'.format(job))

    def command_exe(self, command):
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        return out

    def bowtie2_init(self, fa_file):
        bowtie2 = os.path.join(self.bowtie_root, './bowtie2-build')
        dbname = os.path.join(self.bowtie_root, 'lambda_virus')
        command = '%s -f %s %s'%(bowtie2, fa_file, dbname)
        self.command_exe(command)

    def fa_to_sam(self, fa_file, sam_file):
        print('fasta to sam...')
        bowtie2 = os.path.join(self.bowtie_root, './hisat2')
        dbname = os.path.join(self.bowtie_root, 'fasta')
        # fasta --> -f, fastq --> -q
        command = '%s -x %s -U %s -S %s' % (bowtie2, dbname, fa_file, sam_file)
        self.command_exe(command)
        print('done with fasta to sam')

    def comp_fa_to_sam(self, fa_file, sam_file):
        print('fasta to sam...')
        bowtie2 = os.path.join(self.bowtie_root, './hisat2')
        dbname = os.path.join(self.bowtie_root, 'indexes/genome_tran')

        command = '%s -p 4 -x %s' % (bowtie2, dbname)
        if len(fa_file) < 2:
            return 1

        for i, fa in enumerate(fa_file):
            command = command + ' -{} {}'.format(i + 1, fa)
        command += ' -S {}'.format(sam_file)
        print(command)
        self.command_exe(command)
        print('done with fasta to sam')
        return 0

    def sam_to_bam(self, sam_file):
        print('sam to bam...')
        dirname, fname = os.path.split(sam_file)
        fname = os.path.splitext(fname)[0] + '.bam'
        bam_file = os.path.join(dirname, fname)

        samtool_root = os.path.join(self.root, 'software/samtools-1.9')
        command = '{}/./samtools sort -@ 8 -o {} {}'.format(samtool_root, bam_file, sam_file)
        print(command)
        self.command_exe(command)
        os.remove(sam_file)
        print('done with fasta to bam')

    def auto_run(self, root):
        for path, subdirs, files in os.walk(root):
            for name in files:
                if name.endswith('.bam'):
                    fpath = os.path.join(path, name)
                    self.bam_to_gtf(fpath)
                    self.gff_to_gtf(fpath.replace('.bam', '.gff'))
                    self.gtf_to_db(fpath.replace('.bam', '.gtf'))

    def bam_to_gtf(self, bam_file):
        string_tie_root = os.path.join(self.root, 'software/stringtie-1.3.6')

        print('bam to gff...')
        dirname, fname__ = os.path.split(bam_file)
        fname = os.path.splitext(fname__)[0] + '.gff'
        gtf_file = os.path.join(dirname, fname)

        command = '{str_root}/./stringtie -p 4 -G {str_root}/genes.gtf -o {gtf} -i {bam}' \
                  ''.format(str_root=string_tie_root, gtf=gtf_file, bam=bam_file)
        print(command)
        self.command_exe(command)

        for i in range(2):
            fpath = bam_file.replace('.bam', '_{}.fastq.gz'.format(i + 1))
            if os.path.exists(fpath):
                os.remove(fpath)

        root_dir = '/'.join(bam_file.split('/')[:-1])
        shutil.move(bam_file, os.path.join(root_dir, 'bam', fname__))
        print('done with bam to gtf')

    def bam_to_bed(self, bam_file):
        root = os.path.join(self.root, 'software/samtools')

        print('bam to bed...')
        dirname, fname__ = os.path.split(bam_file)
        fname = os.path.splitext(fname__)[0] + '.bed'
        bed_file = os.path.join(dirname, fname)

        command = '{root}/./convert2bed --input=bam < "{bam}" > "{bed}"'.format(root=root, bam=bam_file, bed=bed_file)
        print(command)
        self.command_exe(command)

        # for i in range(2):
        #     fpath = bam_file.replace('.bam', '_{}.fastq.gz'.format(i + 1))
        #     if os.path.exists(fpath):
        #         os.remove(fpath)

        # root_dir = '/'.join(bam_file.split('/')[:-1])
        # dirname = os.path.join(root_dir, 'bam')
        # if not os.path.exists(dirname):
        #     os.mkdir(dirname)
        # shutil.move(bam_file, os.path.join(dirname, fname__))
        # print('done with bam to gtf')

    def gff_to_gtf(self, gff_file):
        string_tie_root = os.path.join(self.root, 'software/stringtie-1.3.6')

        print('gff to gtf...')
        dirname, fname__ = os.path.split(gff_file)
        fname = os.path.splitext(fname__)[0] + '.gtf'
        gtf_file = os.path.join(dirname, fname)

        command = '{str_root}/./stringtie -p 4 --merge -G {str_root}/genes.gtf -o {gtf} -i {gff}' \
                  ''.format(str_root=string_tie_root, gtf=gtf_file, gff=gff_file)
        print(command)
        self.command_exe(command)

        root_dir = '/'.join(gff_file.split('/')[:-2])
        shutil.move(gff_file, os.path.join(root_dir, 'gff', fname__))
        shutil.move(gtf_file, os.path.join(root_dir, 'gtf', fname))

        print('done with bam to gtf')

    def sub_directory(self, root, ext='.txt'):
        flist = []
        for r, d, f in os.walk(root):
            for file in f:
                if ext in file:
                    flist.append(os.path.join(r, file))
        return flist

    def gtf_to_db(self, root):
        fpaths = self.sub_directory(root, '.gtf')
        db_path = os.path.join(root, 'out', 'RNA_seq.db')
        con = sqlite3.connect(db_path)

        N = len(fpaths)
        for i, fpath in enumerate(fpaths):
            print('{} / {}'.format(i + 1, N))

            dirname, fname = os.path.split(fpath)
            chunksize = 1 << 29  # 512 MB

            if os.path.exists(fpath):
                for df_chunk in pd.read_csv(fpath, sep='\t', chunksize=chunksize, names=self.gtf_columns, comment='#'):
                    df_chunk = df_chunk[df_chunk['feature'] == 'transcript']
                    df_chunk = df_chunk.reset_index(drop=True)
                    df_chunk['chromosome'] = 'chr' + df_chunk['chromosome'].astype('str')
                    attribute = df_chunk['attribute'].str.split('; ')

                    pkm = []
                    for attr in attribute:
                        dict = {}
                        for a in attr[:-1]:
                            key, value = a.split(' ')
                            dict[key] = value.replace('"', '')

                        if 'cov' not in dict or 'gene_name' not in dict:
                            pkm.append([np.nan] * 2)
                            continue
                        pkm.append([dict['gene_name'], dict['cov'].replace(';', '')])

                    df_add = pd.DataFrame(data=pkm, columns=['gene_name', 'cov'])
                    df_con = pd.concat([df_chunk, df_add], axis=1)
                    df_con = df_con.dropna(subset=['gene_name', 'cov'], how='any')
                    tname = os.path.splitext(fname)[0].split('%')[0]
                    df_con.drop(['frame', 'attribute'], axis=1).to_sql(tname, con, index=None, if_exists='append')

    def db_to_bed(self):
        from Database import Database

        dirname = os.path.join(self.root, 'database/Dnase/hg38')
        flist = os.listdir(dirname)
        for fname in flist:
            if fname.endswith('.db'):
                fpath = os.path.join(dirname, fname)
                con = sqlite3.connect(fpath)
                tlist = Database.load_tableList(con)
                for tname in tlist:
                    df = pd.read_sql_query("SELECT * FROM '{}'".format(tname), con)
                    df.to_csv(os.path.join(dirname, tname + '.csv'), sep='\t', index=None, header=False)

    def chain_transfer(self):
        dirname = os.path.join(self.root, 'database/Dnase/hg38')
        flist = os.listdir(dirname)

        for fname in flist:
            if fname.endswith('.csv'):
                root = os.path.join(self.root, 'software/userApps')
                fname, ext = os.path.splitext(fname)
                in_path = os.path.join(dirname, fname + ext)
                out_path = os.path.join(dirname, fname + '_liftover' + ext)
                unmapped = os.path.join(dirname, fname + '_unmap' + ext)

                command = '{root}/./liftOver {in_path} {root}/hg38ToHg19.over.chain.gz {out_path} {unmapped}' \
                          ''.format(root=root, in_path=in_path, out_path=out_path, unmapped=unmapped)
                print(command)
                self.command_exe(command)

    def split_bed_to_db(self):
        dirname = os.path.join(self.root, 'database/Dnase/hg38')
        chunksize = 1 << 20
        fids = {'ENCFF441RET': 'K562', 'ENCFF591XCX': 'HepG2', 'ENCFF716ZOM': 'A549', 'ENCFF775ZJX': 'GM12878',
                'ENCFF912JKA': 'HeLa-S3', 'ENCFF571SSA': 'hESC '}

        flist = [x for x in os.listdir(dirname) if 'liftover' in x]

        for fname__ in flist:
            if fname__.endswith('.csv'):
                fname = os.path.splitext(fname__)[0]
                fid, chromosome, _ = fname.split('_')

                cell_line = fids[fid]
                con = sqlite3.connect(os.path.join(dirname, 'bioinfo_{}.db'.format(cell_line)))
                for df_chunk in pd.read_csv(os.path.join(dirname, fname__), sep='\t', names=['chromosome', 'start', 'end', 'strand'], chunksize=chunksize, low_memory=False):
                    df_chunk = df_chunk[df_chunk['chromosome'].str.len() <= 5]
                    df_chunk.to_sql('{}_{}'.format(fid, chromosome), con, index=None, if_exists='append')

    def bed_to_db(self, dirname):
        chunksize = 1 << 20
        # fids = {'ENCFF441RET': 'K562', 'ENCFF591XCX': 'HepG2', 'ENCFF716ZOM': 'A549', 'ENCFF775ZJX': 'GM12878',
        #         'ENCFF912JKA': 'HeLa-S3', 'ENCFF571SSA': 'hESC', 'ENCFF916WEW': 'MCF7'}

        flist = os.listdir(dirname)
        for fname__ in flist:
            if fname__.endswith('.bed'):
                fname = os.path.splitext(fname__)[0]
                # cell_line = fids[fname]
                # con = sqlite3.connect(os.path.join(dirname, 'bioinfo_{}.db'.format(cell_line)))
                con = sqlite3.connect(os.path.join(dirname, '{}.db'.format(fname)))

                cnt = 0
                for df_chunk in pd.read_csv(os.path.join(dirname, fname__), sep='\t', names=self.bed_columns, chunksize=chunksize, low_memory=False, quotechar=None, quoting=3):
                    print(cnt)
                    df_chunk = df_chunk[['chromosome', 'start', 'end', 'score', 'strand']]
                    cnt += 1

                    df_chr = df_chunk.groupby('chromosome')
                    for chr, df_sub in df_chr:
                        if len(chr) > 5:
                            continue
                        for strand, df_sub_sub in df_sub.groupby('strand'):
                            df_sub_sub.to_sql('{}_hg19_ctss_{}_{}'.format('GM12878', chr, strand), con, index=None, if_exists='append')


if __name__ == '__main__':
    f2b = fa2bed()
    # f2b.bed_to_db()
    # f2b.db_to_bed()
    # f2b.chain_transfer()
    # f2b.split_bed_to_db()

    if f2b.hostname == 'mingyu-Precision-Tower-7810':
        dirname = os.path.join(f2b.root, 'database/Fantom/v5/CAGE/nobacode')
        flist = []
        for r, d, f in os.walk(dirname):
            for file in f:
                if file.endswith('.bam'):
                    f2b.bam_to_bed((os.path.join(r, file)))
        f2b.bed_to_db(dirname)

    else:
        dirname = os.path.join(f2b.root, 'database/Fantom/v5/CAGE')
        flist = []
        for r, d, f in os.walk(dirname):
            for file in f:
                if file.endswith('.bam'):
                    f2b.bam_to_bed((os.path.join(r, file)))
