#! /usr/bin/perl

#####################################################################################################################################################################
### Program:	This program counts if the input nucleotide sequence has any 'N' (case insensitive) characters and the percentage of genome that is 'N'
### Input:	1. Nucleotide sequence (FASTA format) - It can be single sequence or multiple sequences in a single file (multi-fasta). The sequence should be in 
#		continuous single line
###Output:	1. 3 columns as output in the order: Genome size, Number of Ns (case insensitive), percentage of genome with Ns
###Notes:	The program is executed as a job script on UBELIX server using the modules available on the server and can be automated for multiple fasta files
###		hence, run the job batch script named: countNs_job.sh to use this program.
###Author:	Mehta D. (dolly.mehta@unibe.ch)
###Date:	August 17, 2023
######################################################################################################################################################################

#Read the input fasta file
open(IN, $ARGV[0]);
my @lines=<IN>;
close(IN);

#Iterating through each sequence
foreach my $line (@lines)
{
	next if($line=~/^$/);
	
	chomp($line);
	print "$line\t" if($line =~ m/>/);

	#Ignore the header and get the sequence
	if($line !~ m/>/)
	{
		my $n_size=$percentNs=0;
		my $genomeSize=length($line);
		foreach my $char (split(//, $line))
		{
			$n_size+=1 if($char !~ /A|U|T|G|C|a|u|t|g|c/);
		}
		$percentNs=sprintf("%.2f", ($n_size/$genomeSize)*100);
		print "$genomeSize\t$n_size\t$percentNs\n";
	}
}
