#!/bin/bash
#SBATCH --time=4-00:00:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12G
#SBATCH --error=chloro_RNALalifold.err
#SBATCH --output=chloro_RNALalifold.out
#SBATCH --job-name=chloro_RNALalifold
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Note:
#Uses container rnatools_v1.sif located at /storage/workspaces/ips_reislab/reislab/Software/rnatools/. Please make sure to give the path to rnatools_v1.sif in your directory structure

#Commands
date

# Define your file paths
PD=$(pwd)

counter=1;
# Iterate over tasks.
for l in x*; do
	# Calculate the task range for this iteration
	start_task=1
	end_task=$(wc -l < "$i")
	tasks_per_iteration=$end_task;

	i=`echo $l | cut -f2`;
	if [ $counter != 1 ]
	then
		 timer=$(( (counter * iterations) )); 
		 echo -e "Waiting to submit next iteration in $timer seconds\n";
		 sleep $timer;
	fi

	# Submit SLURM array job
	sbatch --array="${start_task}-${end_task}" rnalalifold_array_job.sh "$i";

	((counter++));
done

# Exit after submitting the last iteration to avoid recursion
[ "$SLURM_ARRAY_TASK_ID" == "" ] && exit

date
