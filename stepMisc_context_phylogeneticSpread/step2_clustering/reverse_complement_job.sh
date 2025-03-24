#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=36GB
#SBATCH --error=revComp_plant_chloro.err
#SBATCH --output=revComp_plant_chloro.out
#SBATCH --job-name=revComp_plant_chloro
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

module load vital-it/7
module load UHTS/Analysis/seqtk/1.2

#Commands
in=$1
out=$2

#Get reverse complement of the input sequence - can be multi fasta
seqtk seq -r ${in}> ${out}

#Adding an identifier to distinguish reverse-complementary windows from positive strand windows
sed -i '/^>/s/$/r/' ${out}

date
