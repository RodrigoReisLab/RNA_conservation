#!/bin/bash
#SBATCH --time=4-00:00:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=16GB
#SBATCH --error=mmseq2_chloro_covmode0_PID95_cov80.err
#SBATCH --output=mmseq2_chloro_covmode0_PID95_cov80.out
#SBATCH --job-name=mmseq2_chloro_covmode0_PID95_cov80
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Libraries
module load Anaconda3

#Commands
date
mmseqs easy-cluster all_plant_chloroplast_windows_processed_uniq.fasta mmseq2_covmode0_PID95_cov80 tmp_covmode0_PID95_cov80 -c 0.80 --threads 16 --kmer-per-seq 200 --min-seq-id 0.95 --cov-mode 0 --filter-hits 1

date
