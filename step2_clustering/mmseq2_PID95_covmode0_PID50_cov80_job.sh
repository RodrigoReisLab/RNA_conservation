#!/bin/bash
#SBATCH --time=4-00:00:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=16GB
#SBATCH --error=rna5L_mmseq2_covmode0_PID50_cov80.err
#SBATCH --output=rna5L_mmseq2_covmode0_PID50_cov80.out
#SBATCH --job-name=rna5L_mmseq2_covmode0_PID50_cov80
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Libraries
module load Anaconda3

#Commands
date
mmseqs easy-cluster mmseq2_covmode0_PID95_cov80_rep_seq.fasta mmseq2_repSeqs_covmode0_PID50_cov80 tmp_repSeqs_covmode0_PID50_cov80 -c 0.80 --threads 16 --min-seq-id 0.50 --cov-mode 0 --filter-hits 1

date
