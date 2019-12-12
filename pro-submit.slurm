#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --error=mchajobresults-%J.err
#SBATCH --output=mchajobresults-%J.out
#SBATCH --job-name=mcha_tss_map

# Load modules
echo "Slurm nodes assigned :$SLURM_JOB_NODELIST"
module load R/R-3.6.1-openmpi-4.0.0-ic-2019.3.199
cd /lustre/fs0/home/mcha/Bioinformatics/software/PROmiRNA/external_data

../build/./PROmiRNA -g hg19.fa -c Homo_sapiens.GRCh37.87.gtf -s Hs.five_prime.gff -r hg19_repeats.bed -a hsa.gff3 -m mirna.txt -n mirna_context.txt -p TATA_box_jaspar.psem -w Phastcons/phastCons46way.txt -i bed_files -t 16