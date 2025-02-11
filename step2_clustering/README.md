# This step will help prepare input for RNA structure prediction, remove highly similar windows and group windows based on user-defined sequence similarity thresold.
Perform this step only if you are not sure of the similarity of your input sequences. Incase the sequences are specific (e.g. UTRs of genes of interest), skip this step. You may want to perform a sanity check at the start with a rather quick clustal-omega alignment and making sure the sequences are not identical.
 
**Note:** This step requires use of a conda environments. Kindly use the relevant .yml files provided in this folder.

#### Load conda environment (Anaconda3)
If you are new to conda, you can familiarize yourself with [Manage conda environment] (https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). If you do not have conda installed, please follow [Installation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) guidelines. I have used Anaconda3 for my studies. Once installed, you can activate conda as follows:

	conda activate
	module load Anaconda3 [Specific to users using Ubelix UniBe server]

**Note:** In case the the conda is not initialized in your current shell, then use the following command to set it correctly (Especially on server)

	eval "$(conda shell.bash hook)"

### Optional: Reverse complement sequences using seqtk
The windows (\_uniq.fasta) obtained from Step1 are taken as input. First step is to obtain reverse complementary windows of all the unique windows. This is because the RNA structure prediction step only takes sequences as given and does not automatically handle reverse complementary strand (in genomes- the minus strand). However this step is optional if the input is the desired strand (eg. UTRs of coding-genes will have only a specific strand)

#### Create environment using seqtk_env.yml file.

	conda env create -f seqtk_env.yml

This will create an environment with access to tool seqtk as ready to use. The environment name will be seqtk_env and can be viewed with command:

	conda env list

#### Activate the newly created environment.
	
	conda activate seqtk_env
	
#### Reverse-complement sequences.
i. Use the \_uniq.fasta file obtained from step1 as input 

	seqtk seq -r in_uniq.fasta > out_revComp.fasta

ii. Add identifier to distinguish between forward and reverse complement windows.

	sed -i '/^>/s/$/r/' out_revComp.fasta

iii. Merge the positive and reverse complementary windows in 1 file.

	cat in_uniq.fasta out_revComp.fasta > all_windows.fasta

### Remove highly similar windows and cluster some-what similar windows using mmseqs2
#### Create environment using mmseq2env_list.yml file.

	conda env create -f mmseq2env_list.yml

This will create an environment with access to tool mmseqs2 as ready to use. The environment name will be mmseq2env and can be viewed with command:

	conda env list

#### Activate the newly created environment.
	
	conda activate mmseq2env

#### Remove highly similar sequences

	mmseqs easy-cluster all_plant_chloroplast_windows_processed_uniq.fasta mmseq2_covmode0_PID95_cov80 tmp_covmode0_PID95_cov80 -c 0.80 --threads 16 --kmer-per-seq 200 --min-seq-id 0.95 --cov-mode 0 --filter-hits 1

--min-seq-id 0.95 will help cluster sequences ≥ 95%, which can be modified by the user. Such groups of sequences can be discarded. It also outputs a file ending in \_rep_seq.fasta. This file retains a representative of all these similar sequences.\
Additionally, the other options make sure to obtain appropriate clusters imposing sequence coverage (subject and query) of 80% and --filter-hits 1 is to filter the hits based on the thresholds given.

#### Cluster sequences with some-what sequences similarity
Here, we use the output from the above step (\_rep_seq.fasta) and --min-seq-id is set to 0.50 (50% sequence similarity). Rest all parameters are set as above command.
 
	mmseqs easy-cluster mmseq2_covmode0_PID95_cov80_rep_seq.fasta mmseq2_repSeqs_covmode0_PID50_cov80 tmp_repSeqs_covmode0_PID50_cov80 -c 0.80 --threads 16 --min-seq-id 0.50 --cov-mode 0 --filter-hits 1

Output files will be generated with filenames ending as below, which we will use in the next step:
1. _\_cluster.tsv_: This file includes 2 columns (tab separated). 1st column contains a representative chosen by the tool for a cluster. All members in that cluster are listed in 2nd column. In this pipeline, we use the representative as an identifier for that cluster throughout the study.
Example of the output file:
	
	NC_008590.1_220_38500_38750     NC_031887.1_225_39375_39625
	NC_008590.1_220_38500_38750     NC_044642.1_223_39025_39275
	NC_008590.1_220_38500_38750     NC_029212.1_150_26250_26500r
	NC_008590.1_220_38500_38750     NC_024258.1_224_39200_39450
	NC_008590.1_220_38500_38750     NC_042956.1_218_38150_38400r
	NC_008590.1_220_38500_38750     NC_044800.1_221_38675_38925r
	NC_008590.1_220_38500_38750     NC_027837.1_214_37450_37700r
	NC_008590.1_220_38500_38750     NC_042825.1_93_16275_16525r
	NC_008590.1_220_38500_38750     NC_037832.1_217_37975_38225r
	NC_038124.1_271_47425_47675     NC_038124.1_271_47425_47675
	NC_038124.1_271_47425_47675     NC_036014.1_310_54250_54500r
	NC_038124.1_271_47425_47675     NC_043866.1_305_53375_53625r
	NC_038124.1_271_47425_47675     NC_043822.1_379_66325_66575
	NC_038124.1_271_47425_47675     NC_040966.1_312_54600_54850
	NC_038124.1_271_47425_47675     NC_021647.1_47_8225_8475r


2. _\_all_seqs.fasta_: Contains the sequences in FASTA format separated by an additional header above marked with a cluster representative.

Example output:
	>NC_008590.1_220_38500_38750
	>NC_008590.1_220_38500_38750 
	GAGTAGTAAAGTCTTGTGCTATGAATGCATAAGGAGGTAAAGAGTACATATGTTGAGCTACTAAGGAAGTAATAACCCCTAAAGAAGCTAGAGCAAGACCTAACTGAAAATGAATCGAATTGTTGATTGTGTCATAAAGGCCCTTATGCCCACGCCCTAATCGACCCCCCGGAGGAGTATGCGCTTCTAAAAGATCTTTAATACTGTGCCCAATTCCGAAGTTAGTTCGATACATGTGACCGGCAATGAG
	>NC_042456.1_227_39725_39975 
	AAGTCTTGTGCTATGAATGCATAAGCAGGTAAAGAGTACATGTGTTGAGCTACTAAGGAAGTAATAACCCCTAAAGAAGCTAGAGCAAGGCCTAATTGAAAATGAAGCGAATTATTGATTGTGTCATAAAGACCCTTATGCCCTCGTCCCAATCGTCCCCCCGGAGGAATATGTGCATCTAAAAGATCTTTGATACTGTGTCCAATCCCGAAATTCGTTCGATACATATGACCAGCAACAAGAAAAATAA

#### Preparing cluster member sequences in a parsable manner
Since, the next steps taken multi-fasta file as input and to process each cluster at a time, we split the above \_all_seqs.fasta as each cluster sequence FASTA file with filename as representative to identify the cluster easily.

To identify the number of sequences in a cluster, use the BASH command as follows:
	cut -f1 mmseqs_out_cluster.tsv | sort | uniq -dc

Expected output (tab-separated list): First column contains the number of members, second column contains the identifier (representative of the cluster). We use this to split the cluster fasta file.

	38	NC_008590.1_220_38500_38750
	161	NC_038124.1_271_47425_47675
	186	NC_034645.1_802_140350_140600
	24	NC_041127.1_25_4375_4625
	146	NC_040131.1_286_50050_50300
	64	NC_041473.1_655_114625_114875
	17	NC_028190.1_621_108675_108925
	26	NC_034831.1_51_8925_9175
	36	NC_033910.1_671_117425_117675r
	20	NC_042245.1_489_85575_85825r

For this, the script getClusterSequences.sh helps us obtain the splits and it is executed as follows:

	bash getClusterSequences.sh list

You will get each cluster individual file containing sequences of all its members. Recommended to get this in a separate folder.
