import paramiko
import os


class Server:
    def __init__(self, root, which='newton'):
        self.root = root
        self.which = which
        self.server = '/lustre/fs0/home/mcha/Bioinformatics'

    def connect(self):
        self.ssh = paramiko.SSHClient()
        self.ssh.load_system_host_keys()
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        hostname = '{}.ist.ucf.edu'.format(self.which)
        username = 'mcha'
        password = 'Serenade Seedling cameras Dauphin'
        # key_filename = '/home/mingyu/mcha-keys/mcha_id_rsa_1'
        key_filename = 'D:/mcha-keys/mcha_id_rsa_1'
        self.ssh.connect(hostname, username=username, password=password, key_filename=key_filename)
        print('connected') if self.ssh else print('failed to connect')

    def job_script(self, fname, src_root=None, pversion=3, time='12:00:00'):
        if src_root is None:
            src_root = os.path.join(self.server, 'source/gene_and_mirna')

        script = ['#!/bin/bash', '#SBATCH --nodes=1', '#SBATCH --ntasks-per-node=4', '#SBATCH --cpus-per-task=1',
                  '#SBATCH --time=' + time, '#SBATCH --error=mchajobresults-%J.err',
                  '#SBATCH --output=mchajobresults-%J.out', '#SBATCH --gres=gpu:1',
                  '#SBATCH --job-name=mcha_tss_map\n\n', '# Load modules',
                  'pip install mechanize --user',
                  'echo "Slurm nodes assigned :$SLURM_JOB_NODELIST"',
                  'module load cuda/cuda-9.0',
                  'source {}/venv/bin/activate'.format(src_root),
                  "pip install -r requirements.txt --user",
                  'time python {}'.format('/'.join([src_root, fname]))]
        if self.which == 'stokes':
            script.pop(7)
            script.pop(-3)

        with open('dl-submit.slurm', 'wt', newline='\n') as f:
            f.write('\n'.join(script))

    def upload(self, local_file, server_file):
        print('Upload... ', local_file, ' --> ', server_file)
        tr = self.ssh.get_transport()
        tr.default_max_packet_size = 1 << 32    # 4GB
        tr.default_window_size = 1 << 22        # 4MB

        ftp_client = self.ssh.open_sftp()
        ftp_client.put(local_file, server_file)
        ftp_client.close()
