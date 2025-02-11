#!/bin/bash
#SBATCH --time=4-00:00:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=36GB
#SBATCH --error=drawContext.err
#SBATCH --output=drawContext.out
#SBATCH --job-name=drawContext
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Commands
R CMD BATCH --no-save --no-restore plot_context.R

