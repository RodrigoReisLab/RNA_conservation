# This step is to evaluate the output alignment from step 4 confirming the sequences in the alignment are forming the structure.

This step is divided into multiple sub-steps and it includes the following:  
**i. _trim the alignment_** - alignment from locarna will have the input sequence length (eg. if the input sequences are 100nt - alignment will contain 100nt even if structure is formed of 30nt). Thus, the trimming will help get the relevant 30nt structure with 10-column flank (columns refers to the positions in an alignment) on both sides of the structure.  

**ii.  _evaluate the alignment_** - the trimmed alignment is taken to evaluate using in-house method RNA-ScoRE. Here, each sequence is checked if it can form the structure proposed by LocARNA. It checks for individual base-pairing capability, sequence GC content, sequence redundancy and overall motif formation by the sequence and marks them as 'Passed' or 'Failed' based on these criteria. Number of 'Passed' sequences decide whether the alignment is reliable and are ranked as 'High' (#seqs ≥10), 'Mid' (7 < #seqs > 9) or 'Low' (#seqs < 7).  

**iii. _check for covariation_** - Using R-scape two-tailed test (-s) we check if the alignment marked 'High' or 'Mid' carries any power for covariation. This tool also allows us to check whether the RNAs in the alignment can form conserved alternate structure - if it forms the same structure, the reliability of the motif increases.  

Thus, at the end we will have reliable alignment that we can take forward to the next step.  

### Commands to execute these steps:
**i. _trim the alignment_**

	perl trimAlignment.pl result.stk motif_trimmed.sto motif_trimmed.fasta motif_trimmed.afa

Here, it takes the output alignment from locARNA step (result.stk) and writes 3 output files which can be used for different tools (eg: .sto - Rscape, Infernal; .fasta - RNAz; .afa - SQUARNA). All output files are in alignment format. In case unaligned fasta sequences are required kindly use 'get_unalignedFasta.pl' and 'submit_getUnalignedFasta_job.sh'. The job script can be used as a guide to execute the perl (.pl) program.

**ii. _evaluate the alignment_**

	perl eval_rnaStruct.pl -e _trimmed.sto -d 0 --mt 0.5 -t 0.75 --gc 0.30 motif_trimmed.sto motif_cleaned.sto motif_evaluated.tsv >> allClusters_evaluation.tsv

<pre>
	'e|extension=s'		=> \$ext, 		#any flanking in filename that should be replaced. mandatory to remove .sto from extension - to obtain clean output filename (Here, in this example '_trimmed.sto' from input filename (motif_trimmed.sto) has to be removed.
	'mt|motif_threshold=f'	=> \$motif_thresh, 	#total motif threshold [0.0-1.0] => set to 0.5 for evaluation/cleaning step (column labeled bp_threshold in output file)
	't|bp_threshold=f'	=> \$bp_thresh, 	#each stem in the motif should have x% of base-pairs [0.0-1.0] => set to 0.75 for evaluation/cleaning step. (column labeled HairpinBPs_per in output file)
	'gc|bp_gc=f'		=> \$gc_bp_thresh, 	#each stem in the motif should have x% of GC/CG base-pairs [0.0-1.0] -> set to 0.30
	'd|duplication=i'	=> \$dupl_flag,		#A flag to include duplications in the alignment and all the calculations: 0: duplications OFF - removes duplicationsand 1: duplications ON, retains duplications in the alignment
	'h|help'		=> \$help,
	'manual'		=> \$man
</pre>

It will also output saying if the input alignment is ranked High, Mid or Low on STDOUT. This can be written to any output file (eg: allClusters_evaluation.tsv). If the motif consists of ≤10 base-pairs (ie. total bps ≤ 10) than it is mandatory that the sequence forms all base-pairs. If it does not, the stem-loop is discarded and hence in the evaluation, #BPsActuallyForming will be set to 0.

Passed sequences ranked High and Mid will be written in stockholm alignment format to the output file (eg: 'motif_cleaned.sto')
The tab-separated output file 'motif_evaluated.tsv' will contain information of each sequence in the alignment has passed or failed. The output will contain columns describing the length of the sequence, GC content of the sequence, Threshold deciding how many base-pairs for this motif should be formed, number of total base-pairs, number of opening and closing base-pairs (individual brackets marking the positions as base-paired), percentage of hairpin forming for each stem-loop, number of base-pairs actually formed by the sequence, if there are identical sequences (redundant) and finally if the sequence has passed the evaluation. An output example below shows this motif is a 2 stem-loop structure with 10 base-pairs of which a minimum 5bp should be formed by any sequence.  
Note: TotalBPs, #OpenBPs, #ClosingBPs are calculated from the secondary structure of the alignment (#SS_cons line). Total BPs - 1open and 1 corresponding close bracket - considered as 1basepair. 

| seqID| length | GCcontent | bp_threshold | TotalBPs | #OpenBPs | #ClosingBPs | HairpinBPs_per | #BPsActuallyForming| RedundantRNAcandidates|Evaluation
|:-------------------------------------------:|:------:|:---------:|:------------:|:--------:|:--------:|:-----------:|:--------------:|:--------------------:|:-----------------------------------------:|:--------:|
| NC_035006.1_306_53550_53800r-53986_54172r|94|41|5|10|10|10|60.00,100.00,|5||Passed
| NC_030064.1_279_48825_49075r-49261_49447r|94|33|5|10|10|10|60.00,60.00,|0||Fail
| NC_029813.1_283_49525_49775r-49961_50147r|96|37|5|10|10|10|40.00,60.00,|0||Fail
| NC_041531.1_300_52500_52750r-52936_53122r|96|36|5|10|10|10|60.00,40.00,|0||Fail
| NC_041539.1_294_51450_51700r-51886_52072r|94|43|5|10|10|10|100.00,100.00,|10||Passed
| NC_036977.1_288_50400_50650r-50836_51022r|94|40|5|10|10|10|80.00,100.00,|9||Passed
| NC_034990.1_315_55125_55375r-55561_55747r|94|41|5|10|10|10|80.00,100.00,|9||Passed
| NC_039465.1_65_11375_11625r-11811_11997r|94|39|5|10|10|10|100.00,100.00,|10||Passed
| NC_039924.1_282_49350_49600-49786_49972|104|34|5|10|10|10|80.00,60.00,|4||Fail
| NC_034287.1_275_48125_48375r-48561_48747r|101|42|5|10|8|10|40.00,0.00,|0||Fail
| NC_036052.1_290_50750_51000r-51186_51372r|94|39|5|10|10|10|100.00,100.00,|10||Passed
| NC_040984.1_318_55650_55900r-56086_56272r|94|37|5|10|10|10|80.00,100.00,|9|NC_035566.1_301_52675_52925r-53111_53297r|Passed
| NC_044106.1_304_53200_53450r-53636_53822r|94|39|5|10|10|10|100.00,100.00,|10||Passed
| NC_041266.1_288_50400_50650r-50836_51022r|93|40|5|10|10|10|80.00,80.00,|8||Passed


**iii. _check for covariation_**

	curDir=`pwd`;
	singularity exec -B ${curDir}:/input rnatools_v2.sif bash -c "R-scape -s --cacofold --outdir motif_gc15_twoTest_rscape_out --outname motif_gc15_rscape --voutput --seed 42 /input/motif_cleaned.sto > rscape_motif_twoTest_gc15_out && \
			perl /home/usr/selected_scripts_riboswitch_method/stockholm_to_html.pl /input/motif_cleaned.sto /input/motif_cleaned.html && \
			perl /home/usr/selected_scripts_riboswitch_method/stockholm_to_html.pl /input/motif_gc15_twoTest_rscape_out/motif_gc15_rscape.cacofold.sto /input/motif_gc15_twoTest_rscape_out/motif_gc15_rscape.cacofold.html
		"

R-scape is in the docker container and thus can be used with singularity. _-s_ option performs the two-tailed test and writes the output to a folder 'motif_gc15_twoTest_rscape_out'. This folder should be created before implementing R-scape. To automate this process, use the script 'rscape_eval_twoTest_gc15.sh' as a guide. 

The alignment (motif_cleaned.sto) can be considered as a reliable alignment and used further for the next step.
