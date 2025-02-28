# This folder contains miscellaneous scripts after running the pipeline.

One of the follow-up steps is to see where the RNAs are located in the sequences/ genomes of interest. Thus, a bed file will help with visualizing RNAs in IGV or any gene tracks viewing program.

Thus, the script 'cutsv2bed.pl' will help with converting the file. Here, the motif_evaluated.tsv is taken as input along with corresponding motif_scores.tab and the name of the genome ID (preferably the first column in the motif_evaluated.tsv - either the genome ID or the sequence ID). Note this should not contain co-ordinates. 

e.g. it can be any one of the following:   
NC_000932   
scaffold_1    
scaffold_2    


The script 'cutsv2bed_job.sh' can be used as a guide to execute the perl script. The output BED file (example shown below) is compatible to view in IGV.

<pre>
track name="RNAmotifs" description="RNA motifs identified through extensive computational analyses" visibility=2 colorByStrand="200,90,222 238,187,250" ggfTags=on
NC_000932.1	44865	44913	RNAcoords=NC_000932.1/44865-44913,RNAstrand=minus,SequenceRep=NC_056089.1/44356-44404_-1,model=NC_009259.1_316_55300_55550r,bit_score=50.5,E-value=9.8e-06,Organism=Arabidopsis thaliana 	0	-	44865	44913
NC_000932.1	7786	7837	RNAcoords=NC_000932.1/7786-7837,RNAstrand=plus,SequenceRep=NC_033499.1/7718-7769_1,model=NC_009259.1_316_55300_55550r,bit_score=70.3,E-value=6.3e-10,Organism=Arabidopsis thaliana 	0	+	7786	7837
NC_000932.1	139773	139821	RNAcoords=NC_000932.1/139773-139821,RNAstrand=plus,SequenceRep=NC_066122.1/172374-172422_-1,model=NC_041414.1_478_83650_83900r,bit_score=51.1,E-value=6.3e-06,Organism=Arabidopsis thaliana 	0	+	139773	139821

</pre>
The output BED file is a tab delimited file containing 
	- Column 1 contains Genome ID (This should be identical to the heaer of the genome sequence file loaded in IGV); 
	- Column 2 and 3 containing the RNA co-ordinates 
	- Column 4 containing the description of the RNA. It includes information such as the model from which this motif was identified, the raw bit-score given by cmsearch (how close the homolog is to the model; higher score represents a good homolog) and the e-value (whether the RNA homolog predicted is random, lower the e-value - far from 0, less chance of it being a random hit).
	- Column 5 is supposed to be the score for visualization. Since we do not have any score, we mark it 0
	- Column 6 is the strand where the RNA motif occurs
	- Column 7 and Column 8 are the RNA co-ordinates again (ideally it is used for coding vs non-coding regions like UTRs in the transcripts, which we do not have for RNA motifs).
