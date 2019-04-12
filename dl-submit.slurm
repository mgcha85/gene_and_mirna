#!/bin/bash
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=2
#SBATCH --time=48:00:00
#SBATCH --error=mchajobresults-%J.err
#SBATCH --output=mchajobresults-%J.out
#SBATCH --job-name=mcha_tss_map


# Load modules
echo "Slurm nodes assigned :$SLURM_JOB_NODELIST"
module load anaconda/anaconda3-5.3.0
python Correlation.py