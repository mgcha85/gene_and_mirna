#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --error=mchajobresults-%J.err
#SBATCH --output=mchajobresults-%J.out
#SBATCH --job-name=mcha_tss_map
#SBATCH --gres=gpu:2


# Load modules
echo "Slurm nodes assigned :$SLURM_JOB_NODELIST"
module load cuda/cuda-9.0
module load anaconda/anaconda3
time python /lustre/fs0/home/mcha/Bioinformatics/source/gene_and_mirna/Correlation_gpu.py