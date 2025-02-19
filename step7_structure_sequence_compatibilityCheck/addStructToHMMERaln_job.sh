#! /usr/bin/bash

baseDir=`pwd`;
for i in `cut -f1 listOfRepresentative.txt`;
do
	cd ${baseDir}/${i}_nhmmer;
	echo -e "Adding secondary structure to the hmmer generated alignment for model: ${i}\n";
	perl ${baseDir}/addStructToHMMERaln.pl ${i}_one.sto ${i}_nhmmer.sto > ${i}_nhmmer_withss.sto;
	cd ${baseDir};
done

cd ${baseDir};
