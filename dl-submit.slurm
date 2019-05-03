#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=10
#SBATCH --time=01:30:00
#SBATCH --error=mchajobresults-%J.err
#SBATCH --output=mchajobresults-%J.out
#SBATCH --job-name=mcha_tss_map


# Load modules
echo "Slurm nodes assigned :$SLURM_JOB_NODELIST"
module load anaconda/anaconda3-5.3.0
python /lustre/fs0/home/mcha/Bioinformatics/source/gene_and_mirna/Comparison.py