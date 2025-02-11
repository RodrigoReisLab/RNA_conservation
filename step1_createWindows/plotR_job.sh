#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=36GB
#SBATCH --error=boxplot.err
#SBATCH --output=boxplot.out
#SBATCH --job-name=boxplot
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

module load R/4.1.0-foss-2021a

#Commands
R CMD BATCH --no-save --no-restore boxPlot_distributionOfNsInGenomes.R

date
