#!/bin/bash

IFS=$'\n'
in=$1
echo -e "Chunk considered for getting cluster sequences: $in \n";
while IFS=$'\t' read -r i j
do
	echo -e "Begin to get sequences for cluster: ${j} \n";
	e=$(( 2*i - 1 )); 
	a=$(( $e + 1 ));
	grep -A $e -w -F ">$j" mmseq2_repSeqs_covmode0_PID50_cov80_all_seqs.fasta | tail -n $a > ${j}_cluster.fasta
	echo -e "Cluster sequences for cluster: ${j} written to output file \n";
done<"$in"
