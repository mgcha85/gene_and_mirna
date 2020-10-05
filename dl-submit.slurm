#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --error=mchajobresults-%J.err
#SBATCH --output=mchajobresults-%J.out
#SBATCH --job-name=mcha_tss_map


# Load modules
echo "Slurm nodes assigned :$SLURM_JOB_NODELIST"
source /lustre/fs0/home/mcha/Bioinformatics/source/gene_and_mirna2/venv/bin/activate
time python /lustre/fs0/home/mcha/Bioinformatics/source/gene_and_mirna2/set_go_.py