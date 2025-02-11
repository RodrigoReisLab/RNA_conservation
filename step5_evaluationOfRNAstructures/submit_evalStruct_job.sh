#!/bin/bash
#SBATCH --time=4-00:00:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12G
#SBATCH --error=chloro_evalStruct.err
#SBATCH --output=chloro_evalStruct.out
#SBATCH --job-name=chloro_evalStruct
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Note:
#Uses container rnatools_v1.sif located at /storage/workspaces/ips_reislab/reislab/Software/rnatools/. Please make sure to give the path to rnatools_v1.sif in your directory structure

curDir=`pwd`;
baseDir="/storage/workspaces/ips_reislab/reislab/mehta/chloro_bayer/analyses";
list=$1;

#Checking for any pre-existing output files
if [ -f ${baseDir}/allClusters_evaluation.tsv]
then
	echo -e "Found allClusters_evaluation.tsv in the baseDir. Thus deleting it!\n";
	rm ${baseDir}/allClusters_evaluation.tsv;
fi

for i in `cat $list`;
do
	if [ -f ${baseDir}/${i}/${i}_motif.sto ]
	then
		cd ${curDir}/$i;

		# Evaluating RNAs in the alignment for the motif region (note: the trimmed alignments (output of trimAlignment.pl is used)
		echo -e "Getting the motif region for cluster: $i\n";
		perl ${curDir}/eval_rnaStruct_v1.pl ${i}_motif.sto ${i}_motif_cleaned.sto ${i}_motif_evaluated.tsv >> ${baseDir}/allClusters_evaluation.tsv
		echo -e "Motif region for cluster: $i is evaluated! RNAs that passed the evaluation are written as a cleaned alignment!\n";
	
		cd $curDir;
	else
		echo -e "Check the input file names/ trimmed alignment file and file paths in the code! If no file is found, make sure the input alignment file is in the stockholm format\n";
	fi
done

