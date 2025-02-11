#! usr/bin/bash

#####################################################################################################################################################################
## Program:     Reads fasta multiline input and writes an output fasta in one-line 
## Input:       1. Fasta file (Preferentially in .fa/ .fna) - but not .fasta (since output files are .fasta)
##Output:       1. FASTA one-liner output file (.fasta)

##Notes:        The program is executed as a job script on UBELIX server and hence, 
#               run the job batch script named: fasta_oneLiner_job.sh to use this program.
##Author:       Mehta D. (dolly.mehta@unibe.ch)
##Date:         July 20, 2023
######################################################################################################################################################################

#Read fasta file
fasta_file=$1
basename=$(basename $fasta_file | cut -d. -f1)

if [ ! -e ${basename}_oneLine.fasta ]
then
	while read line ;
	do  
		if [[ ${line:0:1} == ">" ]] 
		then
			echo $'\n' >> ${basename}_oneLine.fasta
			echo $line  >> ${basename}_oneLine.fasta
		else
			echo $line | perl -ne 'chomp and print'>> ${basename}_oneLine.fasta

		fi
	done<${fasta_file}
else
	echo -e "Output file already exists! Delete the file to create a new output file\n";
fi

#Removes an extra empty line at start of the output file!
sed -i '/^$/d' ${basename}_oneLine.fasta
