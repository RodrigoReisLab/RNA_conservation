### This step is to screen for clusters with a potential to form a common structure. 

We do the screening with RNALalifold which is a part of the Vienna RNA package. This package is available within the container (docker or singularity images). Ignore the version of the container used in this step and used the latest version while performing all the steps. Recommended is to have the singularity .sif file in an accessible directory. From this step onwards, we will use this container more frequently.


The cluster FASTA files are taken as input for this step. Most likely in a separate folder named splits. RNALalifold requires inout as an alignment (.aln) file. This can be obtained using mmseqs2 (not tested) or clustal omega. In the studies, we have used clustal omega.

Here, we use 2 input files:
1. a list with clusternames ('cluster_count.tsv' file generated at the end of step2)
2. FASTA sequences of the clusters with paths to its folder.

### To perform clustal omega
Since there are paths to relevant directories, I request to use the clustalo_array_job.sh to execute all the clusters. Please modify the paths to the input data wherever required.
Execute the script as follows:

	bash clustalo_array_job.sh cluster_count.tsv


The output will be .aln file which is used as input for the RNALalifold screening.

**Note:** You will need to have clustal omega tool installed in your system. You can download from [here](http://www.clustal.org/omega/).

### To perform RNALalifold screening
In this step, kindly use the script 'rnalalifold_array_job.sh' as a guide to execute RNALalifold. Within this script you have the following command which uses the container directly.

	singularity exec -B ../${clusterName}:/input /storage/workspaces/ips_reislab/reislab/Software/rnatools/rnatools_v2.1.sif RNALalifold -T 21 --noLP /input/${clusterName}_aligned.aln > ${clusterName}_RNALalifold.out

**Note:** Change path of rnatools_v2.1.sif file to the correct path in your directory.


Options used:
- exec -B: executes a command within the container using folders from the working directory into the container - input directory.
- -T: Temperature setting to fold RNAs.
- --noLP: no lonely pairs allowed (single base-pairs only are not allowed)


Command to execute the script:

	bash rnalalifold_array_job.sh cluster_count.tsv


**Note:** All files and folders are named based on the cluster representative name.

The output will be in the '\_RNALalifold.out' which will include secondary structure in dot-bracket format. Only the clusters with 2 base-pairs predicted are considered to the next step.
Command used to obtain clusters with potential to form secondary structure:

	for i in `find -iname "*_RNALalifold.out"`; do grep -m1 -H "((" $i; done | cut -d':' -f1 | cut -d'/' -f3 | sed -e 's/_RNALalifold.out//g' > RNALalifold_passedList.txt


