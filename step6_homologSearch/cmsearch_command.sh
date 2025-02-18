#!/bin/bash
#SBATCH --time=4-00:00:00
#SBATCH --array=1-620
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=24G
#SBATCH --error=cmsearch_3800genomes_v3_%j.err
#SBATCH --output=cmsearch_3800genomes_v3_%j.out
#SBATCH --job-name=cmsearch_3800genomes_v3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Paths
modelPath="/storage/workspaces/ips_reislab/reislab/mehta/chloro_bayer/v3_afterAnalysesOfLocarnaHits/seedAlns";
PD=$(pwd);
list="high_mid_clustersFrom_v3.txt";
#clusterName=$(awk -v row="${SLURM_ARRAY_TASK_ID}" 'NR==row {print $0}' "$list")

#Commands
for clusterName in `cat $list`
do
	date
	echo -e "Path to models: ${modelPath} \n";
	echo -e "Searching calibrated model: ${clusterName} against NCBI 3800 chloroplast genomes\n";	

	#cmsearch command
	singularity exec -B ${PD}:/input,${modelPath}:/models /storage/workspaces/ips_reislab/reislab/Software/rnatools/rnatools_v2.sif bash -c "
		cmsearch -g -E 0.0001 /models/${clusterName}_motif_cleaned_v3.cm /input/allFasta_oneLine.fa > /input/${clusterName}_hits.txt
		"
	echo -e "Homolog search completed for motif represented by cluster: ${clusterName}\n";
	date
done
