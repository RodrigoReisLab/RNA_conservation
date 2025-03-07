#!/bin/bash
#SBATCH --time=4-00:00:00
#SBATCH --array=1-620
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32G
#SBATCH --error=models620_cmbuild_cmcalibrate.err
#SBATCH --output=models620_cmbuild_cmcalibrate.out
#SBATCH --job-name=models620_cmbuild_cmcalibrate
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Paths to local files and folders
PD=`pwd`;
list="high_mid_clustersFrom_v3.txt"; #Modify the filename accordingly
#clusterName=$(awk -v row="${SLURM_ARRAY_TASK_ID}" 'NR==row {print $0}' "$list")

#Building models from alignments from different cluster
for clusterName in `cat $list`;
do
	singularity exec -B $PD:/input /storage/workspaces/ips_reislab/reislab/Software/rnatools/rnatools_v2.sif bash -c "
		date
		temp=`pwd` && \
		echo -e 'Working in the directory: ${temp} \n' && \
       		cmbuild /input/${clusterName}_motif_cleaned_v3.cm /input/${clusterName}_motif_cleaned_v3.sto && \
       		echo -e 'Covariance model obtained for stockholm(.sto) alignment from v3 evaluated cluster: ${clusterName}\n' && \
       		cmcalibrate /input/${clusterName}_motif_cleaned_v3.cm && \
       		echo -e 'Covariance model ${clusterName} from v3 is calibrated\n' && \
		date
        "
done
