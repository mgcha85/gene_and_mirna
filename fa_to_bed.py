import os
import subprocess
import socket
import pandas as pd


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
        dbname = os.path.join(self.bowtie_root, 'fasta')

        command = '%s -x %s' % (bowtie2, dbname)
        for i, fa in enumerate(fa_file):
            command = command + ' -{} {}'.format(i + 1, fa)
        command += ' -S {}'.format(sam_file)
        self.command_exe(command)
        print('done with fasta to sam')

    def sam_to_bed(self, sam_file, remove=True):
        print('sam to bed...')
        dirname, fname = os.path.split(sam_file)
        fname = os.path.splitext(fname)[0] + '.bed'
        bed_file = os.path.join(dirname, fname)

        samtool_root = os.path.join(self.root, 'software/samtools')
        command = '{}/./sam2bed < {} > {}'.format(samtool_root, sam_file, bed_file)
        self.command_exe(command)
        print('done with fasta to bed')

        columns = ['chromosome', 'start', 'stop', 'id', 'score', 'strand', 'dum1', 'dum2', 'dum3', 'dum4', 'dum5',
                   'dum6', 'dum7', 'dum8', 'dum9', 'dum10', 'flag', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual',
                   'attr']
        df = pd.read_csv(bed_file, sep='\t', names=columns)
        df.loc[:, columns[0]:columns[5]].to_csv(bed_file, sep='\t', index=None, header=False)

        if remove is True:
            os.remove(sam_file)


if __name__ == '__main__':
    f2b = fa2bed()
    fa_file = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/software/bowtie2-2.3.5.1-sra-linux-x86_64/example/reference/SP1.fq'
    sam_file = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/software/bowtie2-2.3.5.1-sra-linux-x86_64/example/reference/SP1.sam'
    f2b.fa_to_sam(fa_file, sam_file)
    f2b.sam_to_bed(sam_file)
