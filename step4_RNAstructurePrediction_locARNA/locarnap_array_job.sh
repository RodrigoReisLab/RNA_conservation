#!/bin/bash
#SBATCH --time=4-00:00:00
#SBATCH --array=1-3136
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12G
#SBATCH --error=chloro_locarnap_%j.err
#SBATCH --output=chloro_locarnap_%j.out
#SBATCH --job-name=chloro_locarnap
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Note:
#Uses container rnatools_v1.sif located at /storage/workspaces/ips_reislab/reislab/Software/rnatools/. Please make sure to give the path to rnatools_v1.sif in your directory structure

#Commands
date
PD=`pwd`;
list=$1;
#clusterName=$(awk -v row="${SLURM_ARRAY_TASK_ID}" 'NR==row {print $1}' "$list")

for clusterName in `cat $list`;
do	
	#Modify the path to the FASTA sequences for each cluster.
	if [ -f /storage/workspaces/ips_reislab/reislab/mehta/chloro_bayer/clustering/splits/${clusterName}_cluster.fasta ]
	then
		mkdir "$clusterName";
		cd "$clusterName";

		echo -e "Input cluster file exists..\nProceeding with further steps.\n";

		#Modify the path to the FASTA sequences for each cluster.
		cp /storage/workspaces/ips_reislab/reislab/mehta/chloro_bayer/clustering/splits/${clusterName}_cluster.fasta .
	
		echo -e "Starting mlocarna run for $clusterName\n";
		singularity exec -B ../${clusterName}:/input /storage/workspaces/ips_reislab/reislab/Software/rnatools/rnatools_v1.sif mlocarna --probabilistic --tgtdir ${clusterName}_locarnap --moreverbose --stockholm --local-progressive --threads=3 --rnafold-temperature=21.0 /input/${clusterName}_cluster.fasta
		echo -e "Finished mlocarna run for $clusterName\n";

		cd "$PD";
	else
		echo -e "$clusterName FASTA file not found in the directory: /storage/workspaces/ips_reislab/reislab/mehta/conservation_project_data/getClusters/mmseq2_clustering/PID50_clusters/splits/ OR \nCheck the input file names and file paths in the code! \n";

	fi
done
date
