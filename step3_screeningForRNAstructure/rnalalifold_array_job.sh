#!/bin/bash
#SBATCH --time=23:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --error=chloro_clusters_rnalalifold_%j.err
#SBATCH --output=chloro_clusters_rnalalifold_%j.out
#SBATCH --job-name=chloro_clusters_rnalalifold
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Note:
#Uses container rnatools_v1.sif located at /storage/workspaces/ips_reislab/reislab/Software/rnatools/. Please make sure to give the path to rnatools_v1.sif in your directory structure

#Commands
date
PD=`pwd`;
list=$1;
clusterName=$(awk -v row="${SLURM_ARRAY_TASK_ID}" 'NR==row {print $2}' "$list")

if [ -f /storage/workspaces/ips_reislab/reislab/mehta/chloro_bayer/RNAlalifold/${clusterName}/${clusterName}_aligned.aln ]
then
	cd ./${clusterName}/
	echo -e "${clusterName} Alignment file present in the folder.\n Proceeding to screen for structures with RNALalifold\n";
	#cp /storage/workspaces/ips_reislab/reislab/mehta/chloro_bayer/RNAlalifold/${clusterName}/${clusterName}_aligned.aln .

	echo -e "Starting RNALalifold run for $clusterName\n";
	singularity exec -B ../${clusterName}:/input /storage/workspaces/ips_reislab/reislab/Software/rnatools/rnatools_v1.sif RNALalifold -T 21 --noLP /input/${clusterName}_aligned.aln > ${clusterName}_RNALalifold.out 
	echo -e "Finished RNALalifold run for $clusterName\n";
	cd $PD
else
	echo -e "$clusterName FASTA file not found in the directory: Check the input file names and file paths in the code! \n";

fi

date
