#! /bin/bash
#SBATCH --time=4-00:00:00
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16GB
#SBATCH --error=cutsv2bed_job.err
#SBATCH --output=cutsv2bed_job.out
#SBATCH --job-name=cutsv2bed_job
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Commands
date
curDir=`pwd`;
if [ ! -d $curDir/context_related ]
then
	mkdir $curDir/context_related;
fi

#cp /storage/workspaces/ips_reislab/reislab/mehta/chloro_bayer/chloroplast_genomes/NCBI_chloroplastID_genomeName.csv .

for i in ` cat full_genomesDB_14034_accessionIDs.txt`;
do
	if [ ! -f $curDir/context_related/${i}_RNAmotifs_gc15.bed ]
	then
		echo -e "$i\n";
		perl cutsv2bed.pl -g $i -s _scores.tab -e _hits_evaluated_gc15.tsv motifsSelected_post_nhmmerStep.txt > $curDir/context_related/${i}_RNAmotifs_gc15.bed
	fi
done
date
