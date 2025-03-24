#!/bin/bash
#SBATCH --time=4-00:00:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32G
#SBATCH --error=chloroHits_rscapeTwotest_gc15.err
#SBATCH --output=chloroHits_rscapeTwotest_gc15.out
#SBATCH --job-name=chloroHits_rscapeTwotest_gc15
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dolly.mehta@unibe.ch

#Note:
#Uses container rnatools_v2.sif located at /storage/workspaces/ips_reislab/reislab/Software/rnatools/. Please make sure to give the path to rnatools_v2.sif in your directory structure
#Note: R-scape is present in the container and can be accessed directly as R-scape

#Commands
date
curDir=`pwd`;
list=$1; #filename passed to the command: homologPassedModels.csv

for i in `cut -f1 $list`
do
	if [ -f ${curDir}/${i}_hits_processed_gc15.sto ] 
	then
		#R-scape
		echo -e "Starting R-scape run with two test for $i alignment file obtained from evaluation with 15% GC threshold\n";
		mkdir ${i}_gc15_twoTest_rscape_out
		singularity exec -B ${curDir}:/input /storage/workspaces/ips_reislab/reislab/Software/rnatools/rnatools_v2.sif bash -c "R-scape -s --cacofold --outdir ${i}_gc15_twoTest_rscape_out --outname ${i}_gc15_rscape --voutput --seed 42 /input/${i}_hits_processed_gc15.sto > rscape_${i}_twoTest_gc15_out && \
			perl /home/usr/selected_scripts_riboswitch_method/stockholm_to_html.pl /input/${i}_hits_processed_gc15.sto /input/${i}_hits_processed_gc15.html && \
			perl /home/usr/selected_scripts_riboswitch_method/stockholm_to_html.pl /input/${i}_gc15_twoTest_rscape_out/${i}_gc15_rscape.cacofold.sto /input/${i}_gc15_twoTest_rscape_out/${i}_gc15_rscape.cacofold.html
		"
		echo -e "Finished R-scape run two test for $i for alignment files obtained from evaluation with 15% GC threshold\n";

		cd "$curDir";
	else
		echo -e "Check the input file names and file paths in the code! \n";
	fi
done
date
