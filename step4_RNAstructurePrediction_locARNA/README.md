# This step will find a common structure in the sequences given as input. We use a progressive aligner - LocARNA.

Here we use the list of clusters which passed RNALalifold step and find what the common structure that the sequences share. We perform this using LocARNA which performs sequence and structure guided predictions. This is the crucial and the first step to obtain the conserved structure between sequences.

The script 'locarnap_array_job.sh' can be used as a guide to perform LocARNA which uses unaligned sequences with paths to the correct input files and folders. Thus, use the same input FASTA sequences generated in Step2. LocARNA is in the container - the version is only used here as an example and kindly use the latest version. 

Command to execute locarna from within container. Modify the path to the container sif file. Modify the paths in the script file as well.

	singularity exec -B ../${clusterName}:/input /storage/workspaces/ips_reislab/reislab/Software/rnatools/rnatools_v2.1.sif mlocarna --probabilistic --tgtdir ${clusterName}_locarnap --moreverbose --stockholm --local-progressive --threads=3 --rnafold-temperature=21.0 /input/${clusterName}_cluster.fasta  

Options used:
- --probabilistic: the alignments are scored based on match probabilities 
- --tgtdir: the name of the output directory created new
- --moreverbose: writes all the details of the program executed in a log file. It is mainly to help troubleshoot, if any problems.
- --stockholm: writes the output alignment in stockholm format. It is a specific format with sequence aligned and mapped with base-pairing positions in WUSS format.
- --local-progressive: performs local alignment which aids in finding local conserved RNA structures rather than long-range interactions within 250nt windows. Progressive alignment helps adding sequences one by one to the alignment.
- --rnafold-temperature: temperature for folding each sequence can be set as per the user.

The script can be executed as

	bash locarnap_array_job.sh RNALalifold_passedList.txt


This will generate an output folder named '\_locarnap'. Within this folder, there will be alignment file - 'result.stk'. This is the conserved RNA structure predicted for the sequences in the cluster.



