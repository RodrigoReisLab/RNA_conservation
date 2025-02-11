#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12G
#SBATCH --error=chloro_clustalo_%j.err
#SBATCH --output=chloro_clustalo_%j.out
#SBATCH --job-name=chloro_clustalo
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Note:
#Uses container rnatools_v1.sif located at /storage/workspaces/ips_reislab/reislab/Software/rnatools/. Please make sure to give the path to rnatools_v1.sif in your directory structure

#Modules
module load vital-it/7
module load SequenceAnalysis/MultipleSequenceAlignment/clustal-omega/1.2.4

#Commands
date
PD=`pwd`;
list=$1;

echo -e "Considering chunk: $list\n" > clustalo_${list}.log;
for i in `cut -f2 "$list"`;
do
	clusterName=$i;
	if [ -f /storage/workspaces/ips_reislab/reislab/mehta/chloro_bayer/clustering/splits/${clusterName}_cluster.fasta ]
	then
		mkdir "$clusterName";
		cd "$clusterName";

		echo -e "Input cluster file exists..\nProceeding with further steps.\n" >> clustalo_${list}.log;
		cp /storage/workspaces/ips_reislab/reislab/mehta/chloro_bayer/clustering/splits/${clusterName}_cluster.fasta .

		echo -e "Starting clustalo run for $clusterName\n" >> clustalo_${list}.log;
		clustalo -i ${clusterName}_cluster.fasta --percent-id --full --distmat-out ${clusterName}_distMat.csv -o ${clusterName}_aligned.aln --outfmt clu
		echo -e "Finished clustalo run for $clusterName\n" >> clustalo_${list}.log;

		cd "$PD";
	else
		echo -e "$clusterName FASTA file not found in the directory: /storage/workspaces/ips_reislab/reislab/mehta/chloro_bayer/clustering/splits/ OR \nCheck the input file names and file paths in the code! \n" >> clustalo_${list}.log;

	fi
done

date
