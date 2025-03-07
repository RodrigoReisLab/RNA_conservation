#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=24GB
#SBATCH --error=windows.err
#SBATCH --output=windows.out
#SBATCH --job-name=createOverlappingWindows
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Commands
for i in `ls *_oneLine.fasta | cut -d'.' -f1`
do 
	perl createWindows.pl -f ${i}.fasta -w 250 -p 75 -O ${i}_windows.fa
done

date
