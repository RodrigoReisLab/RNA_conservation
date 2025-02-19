#Performed locally with recommended R-scape version: rscape_v2.0.4.a
#Perform the following command to get the fasta unaligned sequences prior to running this script:
# ~/Documents/rscape_v2.0.4.a/lib/hmmer/easel/miniapps/esl-reformat -o ${clusterName}.fa fasta ${clusterName}.sto

#Commands
date
PD=`pwd`;
list="listOfRepresentative.txt";
#clusterName=$(awk -v row="${SLURM_ARRAY_TASK_ID}" 'NR==row {print $1}' "$list")

for i in `cut -f1 $list`;
do
	clusterName=$i;
	if [ -f ${PD}/${clusterName}.sto ]
	then
	        if [[ ! -d ${PD}/${clusterName}_nhmmer ]]
	        then
	                mkdir ${PD}/${clusterName}_nhmmer;
	        fi
	
	        cd ${PD}/${clusterName}_nhmmer
	
	        #Getting required information and files
	        grep "$clusterName" ${PD}/listOfRepresentative.txt | cut -f2 > ${PD}/${clusterName}_nhmmer/list.txt
	        rep=`grep "$clusterName" ${PD}/listOfRepresentative.txt | cut -f2`
	        cp ${PD}/${clusterName}.sto .
	        cp ${PD}/${clusterName}.fa .

	        #Reformatting Infernal ouput
	        echo -e "Performing nhmmer and related scripts for model: ${clusterName}\n";
	        ~/Documents/rscape_v2.0.4.a/lib/hmmer/easel/miniapps/esl-sfetch ${clusterName}.sto $rep > ${clusterName}_one.fa
	        ~/Documents/rscape_v2.0.4.a/lib/hmmer/easel/miniapps/esl-alimanip --seq-k list.txt ${clusterName}.sto > ${clusterName}_one.v0.sto
	        ~/Documents/rscape_v2.0.4.a/lib/hmmer/easel/miniapps/esl-alimask --gapthresh 0.9 -g ${clusterName}_one.v0.sto > ${clusterName}_one.sto
	        echo -e "Input files prepared for performing nhmmer for cluster: ${clusterName}\nPerforming nhmmer..!!\n"
	        ~/Documents/rscape_v2.0.4.a/lib/hmmer/src/hmmbuild ${clusterName}_one.hmm ${clusterName}_one.fa
	        ~/Documents/rscape_v2.0.4.a/lib/hmmer/src/hmmsearch -A ${clusterName}_nhmmer.sto ${clusterName}_one.hmm ${clusterName}.fa > ${clusterName}_hmmsearch.log
	       	echo -e "nhmmer step completed!\nRunning R-scape with cacofold..\n";
		if [[ ! -d  ${clusterName}_rscape_nhmmer ]]
		then
		        mkdir ${clusterName}_rscape_nhmmer
		fi
	        ~/Documents/rscape_v2.0.4.a/bin/R-scape --cacofold --outdir ${clusterName}_rscape_nhmmer ${clusterName}_nhmmer.sto > ${clusterName}_rscape.log
	        echo -e "nhmmer and Rscape completed for ${clusterName}\n";
	        date

fi

cd $PD;
done
