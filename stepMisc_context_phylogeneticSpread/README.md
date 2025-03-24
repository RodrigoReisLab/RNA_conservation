# This folder contains miscellaneous scripts after running the pipeline.

One of the follow-up steps is to see where the RNAs are located in the sequences/ genomes of interest. Thus, a bed file will help with visualizing RNAs in IGV or any gene tracks viewing program.

Thus, the script 'cutsv2bed.pl' will help with converting the file. Here, the motif_evaluated.tsv is taken as input along with corresponding motif_scores.tab and the name of the genome ID (preferably the first column in the motif_evaluated.tsv - either the genome ID or the sequence ID). Note this should not contain co-ordinates. 

e.g. it can be any one of the following:
NC_000932
scaffold_1  
scaffold_2  


The script 'cutsv2bed_job.sh' can be used as a guide to execute the perl script. 

Add about hot to obtain for IGV
