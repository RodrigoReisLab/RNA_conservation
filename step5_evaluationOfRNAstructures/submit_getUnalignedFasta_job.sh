#!/bin/bash
#SBATCH --time=4-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12G
#SBATCH --error=getUnalignedFasta.err
#SBATCH --output=getUnalignedFasta.out
#SBATCH --job-name=getUnalignedFasta
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Note:
#Uses container rnatools_v1.sif located at /storage/workspaces/ips_reislab/reislab/Software/rnatools/. Please make sure to give the path to rnatools_v1.sif in your directory structure

curDir=`pwd`;
list=$1;

for i in `cat $list`;
do
	if [ -f ${curDir}/${i}/${i}_motif.fasta ]
	then

		cd ${curDir}/${i}/
		# Removing gaps from the FASTA files
		echo -e "Found FASTA file for cluster: $i\n";
		perl ${curDir}/get_unalignedFasta.pl ${i}_motif.fasta ${i}_motif_ungapped.fa

		cd $curDir;
	else
		echo -e "Could not find the FASTA file for the cluster: $i\n Hint: Check the output of trimAlignment.pl for the cluster  OR \nCheck the input file names and file paths in the code! \n";
	fi
done

