#!/bin/bash
#SBATCH --time=4-00:00:00
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=64GB
#SBATCH --error=bayer_fasta_oneLine.err
#SBATCH --output=bayer_fasta_oneLine.out
#SBATCH --job-name=bayer_convertInputToOneLineFasta
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Commands
for i in `ls *.fa`
do 
	bash fasta_oneLiner.sh $i
done

date
