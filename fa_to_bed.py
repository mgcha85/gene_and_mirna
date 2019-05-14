import os
import subprocess
import socket
import pandas as pd
import sqlite3


class fa2bed:
    def __init__(self):
        hostname = socket.gethostname()
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
            # self.root = 'D:/Bioinformatics'
        elif hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'
        self.bowtie_root = os.path.join(self.root, 'software/hisat2-2.1.0')
        self.gtf_columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

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

        command = '%s -x %s' % (bowtie2, dbname)
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
                    self.gtf_to_db(fpath.replace('.bam', '.gtf'))

    def bam_to_gtf(self, bam_file):
        string_tie_root = os.path.join(self.root, 'software/stringtie-1.3.3b')

        print('bam to gtf...')
        dirname, fname = os.path.split(bam_file)
        fname = os.path.splitext(fname)[0] + '.gtf'
        gtf_file = os.path.join(dirname, fname)

        command = '{str_root}/./stringtie -p 8 -G {str_root}/genes.gtf -o {gtf} {bam}'.format(str_root=string_tie_root,
                                                                                        gtf=gtf_file, bam=bam_file)
        print(command)
        self.command_exe(command)

        os.remove(bam_file)
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
                        for a in attr:
                            key, value = a.split(' ')
                            dict[key] = value.replace('"', '')
                        pkm.append([dict['FPKM'], dict['TPM'].replace(';', '')])

                    df_add = pd.DataFrame(data=pkm, columns=['FPKM', 'TPM'])
                    df_con = pd.concat([df_chunk, df_add], axis=1)
                    tname = os.path.splitext(fname)[0].split('%')[0]
                    df_con.drop(['frame', 'attribute'], axis=1).to_sql(tname, con, index=None, if_exists='append')

    def sam_to_bed(self, sam_file, to_db=True, remove=True):
        print('sam to bed...')
        dirname, fname = os.path.split(sam_file)
        fname = os.path.splitext(fname)[0] + '.bed'
        bed_file = os.path.join(dirname, fname)

        # samtool_root = os.path.join(self.root, 'software/samtools')
        # command = '{}/./convert2bed --input=SAM --output=BED < {} > {}'.format(samtool_root, sam_file, bed_file)
        # self.command_exe(command)
        # print('done with fasta to bed')

        self.gtf_to_db(gtf_file)

        if remove is True:
            os.remove(sam_file)


if __name__ == '__main__':
    f2b = fa2bed()
    f2b.gtf_to_db('/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/database/RNA-seq')
