# This step is to evaluate the output alignment from step 4 confirming the sequences in the alignment are forming the structure.

This step is divided into multiple sub-steps and it includes the following:
i. _trim the alignment_ - alignment from locarna will have the input sequence length (eg. if the input sequences are 100nt - alignment will contain 100nt even if structure is formed of 30nt). Thus, the trimming will help get the relevant 30nt structure with 10-column flank (columns refers to the positions in an alignment) on both sides of the structure.
ii.  _evaluate the alignment_ - the trimmed alignment is taken to evaluate using in-house method RNA-ScoRE. Here, each sequence is checked if it can form the structure proposed by LocARNA. It checks for individual base-pairing capability, sequence GC content, sequence redundancy and overall motif formation by the sequence and marks them as 'Passed' or 'Failed' based on these criteria. Number of 'Passed' sequences decide whether the alignment is reliable and are ranked as 'High' (#seqs ≥10), 'Mid' (7 < #seqs > 9) or 'Low' (#seqs < 7).
iii. _check for covariation_ - Using R-scape two-tailed test (-s) we check if the alignment marked 'High' or 'Mid' carries any power for covariation. This tool also allows us to check whether the RNAs in the alignment can form conserved alternate structure - if it forms the same structure, the reliability of the motif increases.

Thus, at the end we will have reliable alignment that we can take forward to the next step. 

### Commands to execute these steps:
i. _trim the alignment_

	perl trimAlignment.pl result.stk ${i}_motif.sto ${i}_motif.fasta ${i}_motif.afa

Here, it takes the output alignment from locARNA step (result.stk) and writes 3 output files which can be used for different tools (eg: .sto - Rscape, Infernal; .fasta - RNAz; .afa - SQUARNA).  

ii. _evaluate the alignment_

	perl eval_rnaStruct_v3.pl -e _motif.sto -d 0 --mt 0.5 -t 0.75 --gc 0.30 motif_trimmed.sto motif_cleaned.sto motif_evaluated.tsv >> allClusters_evaluation.tsv


