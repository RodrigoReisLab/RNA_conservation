#!/bin/bash
#SBATCH --time=4-00:00:00
#SBATCH --array=1-3
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=16GB
#SBATCH --error=getClusterSeqs_PID50_%j.err
#SBATCH --output=getClusterSeqs_PID50_%j.out
#SBATCH --job-name=getClusterSeqs_PID50_arr
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Commands
date
echo -e "Starting job..!!\n";
cd /storage/workspaces/ips_reislab/reislab/mehta/chloro_bayer/clustering/splits
config=./config.txt
chunk=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

echo -e "Processing $chunk\n";
bash getClusterSequences.sh $chunk
echo -e "Got all the clusters sharing 50% identity..!!\nJob has ended..!!\n";

date
