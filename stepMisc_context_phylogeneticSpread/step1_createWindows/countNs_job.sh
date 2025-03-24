#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=36GB
#SBATCH --error=countNs.err
#SBATCH --output=countNs.out
#SBATCH --job-name=countNs
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Commands
for i in `ls *_oneLine.fasta | cut -d'.' -f1`
do
	if [ ! -s ${i}_nuclInfo.csv ]
	then
		perl countNs.pl ${i}.fasta > ${i}_nuclInfo.csv
	fi
done


date
