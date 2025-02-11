#!/bin/bash
#SBATCH --time=4-00:00:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=12G
#SBATCH --error=chloro_trimAln.err
#SBATCH --output=chloro_trimAln.out
#SBATCH --job-name=chloro_trimAln
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Note:
#Uses container rnatools_v1.sif located at /storage/workspaces/ips_reislab/reislab/Software/rnatools/. Please make sure to give the path to rnatools_v1.sif in your directory structure

#Modules (Comment/Uncomment accordingly)
#On submit02 node
module load HMMER/3.3.2-gompi-2021a

#On submit021/03 node
#module load vital-it/7
#module load SequenceAnalysis/HMM-Profile/hmmer/3.1b2

curDir=`pwd`;
baseDir="/storage/workspaces/ips_reislab/reislab/mehta/chloro_bayer/locarna";
list=$1;

for i in `cat $list`;
do
	if [ -f ${baseDir}/${i}/${i}_locarnap/results/result.stk ]
	then
		mkdir ${curDir}/$i;
		cd ${curDir}/$i;

		# Copying the relevant alignment file (Output by locarna)
		cp ${baseDir}/${i}/${i}_locarnap/results/result.stk ${curDir}/$i;

		# Trimming to get the motif region from the entire alignment
		echo -e "Getting the motif region for cluster: $i\n";
		perl ${curDir}/trimAlignment.pl result.stk ${i}_motif.sto ${i}_motif.fasta ${i}_motif.afa
		esl-reformat --informat stockholm -o ${i}_motif.aln clustal ${i}_motif.sto
		echo -e "Motif region for cluster: $i written as stockholm, fasta and clustalw format\n";
	
		cd $curDir;
	else
		echo -e "Locarna has not found any conserved structures among sequenes in the cluster: $i  OR \nCheck the input file names and file paths in the code! \n";
	fi
done

