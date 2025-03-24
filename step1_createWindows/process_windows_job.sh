#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=512GB
#SBATCH --error=removeNs_polyN_windows.err
#SBATCH --output=removeNs_polyN_windows.out
#SBATCH --job-name=processingWindows
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Commands
echo -e "Starting job: Removal of sequences with Ns and poly N sequences!\n";
perl removeNs_polyN_windows.pl plant_chloroplast_windows.fasta plant_chloroplast_windows_processed
echo -e "Removed sequences with Ns and poly N asequences!\n";

echo -e "Starting job: Removal of duplicate sequences!\n";
perl removeDuplicates.pl plant_chloroplast_windows_processed_noNs_polyN
echo -e "Removed duplicated sequences and unique sequences are written to output file!\n";

date
