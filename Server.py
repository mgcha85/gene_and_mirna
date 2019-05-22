import paramiko
import socket
import os
from download_RNA_seq import Download_RNA_seq
from fa_to_bed import fa2bed
import time


class Server:
    def __init__(self):
        hostname = socket.gethostname()
        if hostname == 'mingyu-Precision-Tower-7810':
            self.root = '/media/mingyu/70d1e04c-943d-4a45-bff0-f95f62408599/Bioinformatics'
        elif hostname == 'DESKTOP-DLOOJR6':
            self.root = 'D:/Bioinformatics'
        elif hostname == 'mingyu-Inspiron-7559':
            self.root = '/media/mingyu/8AB4D7C8B4D7B4C3/Bioinformatics'
        else:
            self.root = '/lustre/fs0/home/mcha/Bioinformatics'

        self.server = '/lustre/fs0/home/mcha/Bioinformatics'
        self.rna_dir = os.path.join(self.server, 'database/RNA-seq/fastq/{}')
        self.f2b = fa2bed()

    def connect(self, which='newton'):
        self.ssh = paramiko.SSHClient()
        self.ssh.load_system_host_keys()

        hostname = '{}.ist.ucf.edu'.format(which)
        username = 'mcha'
        password = 'Serenade Seedling cameras Dauphin'
        key_filename = '/home/mingyu/mcha-keys/mcha_id_rsa_1'
        self.ssh.connect(hostname, username=username, password=password, key_filename=key_filename)

    def job_script(self, fname, time='12:00:00', which='newton'):
        src_root = os.path.join(self.server, 'source/gene_and_mirna')
        script = ['#!/bin/bash', '#SBATCH --nodes=4', '#SBATCH --ntasks-per-node=4', '#SBATCH --time='+time,
                  '#SBATCH --error=mchajobresults-%J.err', '#SBATCH --output=mchajobresults-%J.out',
                  '#SBATCH --gres=gpu:1','#SBATCH --job-name=mcha_tss_map\n\n', '# Load modules',
                  'echo "Slurm nodes assigned :$SLURM_JOB_NODELIST"',
                  'module load cuda/cuda-9.0', 'module load anaconda/anaconda3',
                  'time python {}'.format(os.path.join(src_root, fname))]
        # if which == 'stokes':
        script.pop(6)
        script.pop(-3)

        with open('dl-submit.slurm', 'wt') as f:
            f.write('\n'.join(script))

    def upload(self, local_file, server_file):
        print('Upload... ', local_file, ' --> ', server_file)
        tr = self.ssh.get_transport()
        tr.default_max_packet_size = 1 << 32    # 4GB
        tr.default_window_size = 1 << 22        # 4MB

        ftp_client = self.ssh.open_sftp()
        ftp_client.put(local_file, server_file)
        ftp_client.close()

    def run(self):
        for i in range(1, 2):
            if i < 5:
                which = 'stokes'
            else:
                which = 'newton'
            self.connect(which)
            self.job_script('download_RNA_seq.py', which)

            src_root = os.path.join(self.server, 'source/gene_and_mirna')
            self.upload('dl-submit.slurm', os.path.join(src_root, 'dl-submit.slurm'))
            self.upload('requirements.txt', os.path.join(src_root, 'requirements.txt'))

            dirname = self.rna_dir.format(i)
            print(dirname)
            stdin, stdout, stderr = self.ssh.exec_command("cd {};sbatch {}/dl-submit.slurm".format(dirname, src_root))
            job = stdout.readlines()[0].replace('\n', '').split(' ')[-1]
            print('job ID: {}'.format(job))


if __name__ == '__main__':
    server = Server()
    server.run()
