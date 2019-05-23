import pandas as pd
import os
import subprocess
import socket
from Server import Server
import sys
import shutil


class Convert:
    def __init__(self):
        self.hostname = socket.gethostname()
        if self.hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif self.hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif self.hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

    def to_server(self):
        server = Server()
        server.connect()

        local_path = sys.argv[0]
        dirname, fname = os.path.split(local_path)

        server.job_script(fname, time='03:00:00')

        server_root = os.path.join(server.server, 'source/gene_and_mirna')
        server_path = local_path.replace(dirname, server_root)

        server.upload(local_path, server_path)
        server.upload('dl-submit.slurm', os.path.join(server_root, 'dl-submit.slurm'))

        stdin, stdout, stderr = server.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(server_root, server_root))
        job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
        print('job ID: {}'.format(job))

    def command_exe(self, command):
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        return out

    def bam_to_gtf(self, bam_file):
        string_tie_root = os.path.join(self.root, 'software/stringtie-1.3.6')

        print('bam to gff...')
        dirname, fname = os.path.split(bam_file)
        fname = os.path.splitext(fname)[0] + '.gff'
        gff_file = os.path.join(dirname, fname)

        command = '{str_root}/./stringtie -p 8 -G {str_root}/genes.gtf -o {gff} -i {bam}' \
                  ''.format(str_root=string_tie_root, gff=gff_file, bam=bam_file)
        print(command)
        self.command_exe(command)

        print('done with bam to gtf')

    def gff_to_gtf(self, gff_file):
        string_tie_root = os.path.join(self.root, 'software/stringtie-1.3.6')

        print('gff to gtf...')
        dirname, fname = os.path.split(gff_file)
        fname = os.path.splitext(fname)[0] + '.gtf'
        gtf_file = os.path.join(dirname, fname)

        command = '{str_root}/./stringtie --merge -e -p 8 -G {str_root}/genes.gtf -o {gtf} -i {gff}' \
                  ''.format(str_root=string_tie_root, gtf=gtf_file, gff=gff_file)
        print(command)
        self.command_exe(command)
        print('done with bam to gtf')

    def arrange_files(self, bam_file):
        gff_file = bam_file.replace('.bam', '.gff')
        gtf_file = bam_file.replace('.bam', '.gtf')

        dirname, fname__ = os.path.split(bam_file)
        fname, ext = os.path.splitext(fname__)
        super_dir = '/'.join(dirname.split('/')[:-2])
        
        for from_path, type in zip([bam_file, gff_file, gtf_file], ['bam', 'gff', 'gtf']):
            to_path = os.path.join(super_dir, type, fname + '.' + type)
            shutil.move(from_path, to_path)

    def split_files(self, root_dir, batch_size, ext='.fastq'):
        flist = os.listdir(root_dir)
        flist = [os.path.join(root_dir, x) for x in flist if x.endswith(ext)]

        N = len(flist)
        M = (N + (batch_size - 1)) // batch_size
        for i, fpath in enumerate(flist):
            dirnum = (i // M) + 1
            _, fname = os.path.split(fpath)

            dirname = os.path.join(root_dir, str(dirnum))
            if not os.path.exists(dirname):
                os.mkdir(dirname)
            shutil.move(fpath, os.path.join(dirname, fname))

    def run(self):
        # dirname = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics/database/RNA-seq/test'
        dirname = os.getcwd()

        flist = os.listdir(dirname)
        fpaths = [os.path.join(dirname, x) for x in sorted(flist) if x.endswith('.bam')]
        for fpath in fpaths:
            self.bam_to_gtf(fpath)
            self.gff_to_gtf(fpath.replace('.bam', '.gff'))
            self.arrange_files(fpath)


if __name__ == '__main__':
    con = Convert()
    con.run()
