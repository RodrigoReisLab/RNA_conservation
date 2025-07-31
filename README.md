# **Genome-wide discovery of conserved RNA secondary structures**
-----------------------------------------------------------------

_Copyright 2025 [Dolly Mehta](https://github.com/DollyCM)_

Please cite: [preprint_soon](https://www.biorxiv.org/)

This pipeline was created to find conserved RNA structures in genomes or sequence datasets (eg. UTR regions, RNAs from transcriptomes). 

We have applied this pipeline to [plastid genomes](https://www.biorxiv.org/).

## **Pipeline to predict conserved RNA structures in genomes or sequences of interest.**

The pipeline is in 6-7 steps starting from input of sequences. Please make sure the sequences are in the DNA format (no 'U'), if starting from step 0. 

We recommend the use of high performance computing cluster (HPC). Our scripts were made for a SLURM-managed HPC, but they can easily be adaptated to others.

<div align="center">
<img src="pipeline.png", width="500px">
</div>
Figure 1: Overview of the pipeline. Each of the box are color-coded and marked to help identify the steps to be followed. Each of the steps are dispensable depending on which step user wishes to perform. Scripts and tools relevant to each step are in the respective directories.


### For the entire pipeline, following containers and environments files will be easy to use.

1. rnatools Docker file -  easy access to tools like Vienna RNA package, LocARNA, R-scape and HMMER. It can be taken from [dockerhub.](https://hub.docker.com/r/dollycm/rnatools) and use the following commands to get the container. 
 You can use either docker or singularity [Recommended] to pull the latest image. The current version is v2.1.

        docker pull dollycm/rnatools:v2.1
        
     OR

        singularity pull docker://dollycm/rnatools:v2.1

2. mmseqs2 YAML - for clustering of sequences. YAML file is provided in the respective Step2 clustering folder

3. seqtk YAML - for basic sequence manipulations, including reverse complement input sequences. YAML file is provided in the respective Step2 clustering folder.

4. RNA-SCoRE program - for evaluating alignments for structure. It can be obtained from [here](https://github.com/RodrigoReisLab/RNA-SCoRE)

StepMisc folder contains miscellaneous steps to use the motifs for analyses in obtaining phylogenetic spread and gene context.

**Note:** Some of the job bash scripts have commands for SLURM based execution. If executing the scripts on local machine, directly run the programs as normal bash scripts. The SLURM commands will not affect the performance of the scripts.
