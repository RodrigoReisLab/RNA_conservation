#!/bin/bash
#SBATCH --time=4-00:00:00
#SBATCH --array=1-620
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --error=reformat_cmalign_%j.err
#SBATCH --output=reformat_cmalign_%j.out
#SBATCH --job-name=reformat_cmalign
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Note:
#Uses container rnatools_v2.sif located at /storage/workspaces/ips_reislab/reislab/Software/rnatools/. Please make sure to give the path to rnatools_v2.sif in your directory structure. The sif file is generated from docker image from dockerhub: https://hub.docker.com/repository/docker/dollycm/rnatools/general

#Commands
date
PD=`pwd`;
modelPath="/storage/workspaces/ips_reislab/reislab/mehta/chloro_bayer/v3_afterAnalysesOfLocarnaHits/seedAlns"
list="high_mid_clustersFrom_v3.txt";
#clusterName=$(awk -v row="${SLURM_ARRAY_TASK_ID}" 'NR==row {print $0}' "$list")

for clusterName in `cat $list`;
do
	if [ -f ${PD}/${clusterName}_hits.txt ]
	then
		#Reformatting Infernal ouput
		echo -e "Processing Infernal output to get parseable findings for ${clusterName}\n";
		singularity exec -B $PD:/input,$modelPath:/models /storage/workspaces/ips_reislab/reislab/Software/rnatools/rnatools_v2.sif bash -c "
	      	perl /home/usr/selected_scripts_riboswitch_method/cmsearch_reformatv1_1.pl -s ${clusterName}_scores.tab /input/${clusterName}_hits.txt /input/${clusterName}_hits.fna && \
	       	echo -e 'Using cmalign to align identified hits to the cluster model: ${clusterName}\n' && \
	       	cmalign -o /input/${clusterName}_hits.sto /models/${clusterName}_motif_cleaned_v3.cm /input/${clusterName}_hits.fna && \
	       	echo -e 'Writing an html output for stockholm(.sto) alignment representing hits from cluster: ${clusterName}\n' && \
	       	perl /home/usr/selected_scripts_riboswitch_method/stockholm_to_html.pl /input/${clusterName}_hits.sto /input/${clusterName}_hits.html
	        "
		echo -e "Hits identified from Infernal are processed and output as scores (.tab), fasta (.fna), stockholm (.sto) and colored alignment view (.html) for ${clusterName}\n";
		date
	fi
done

