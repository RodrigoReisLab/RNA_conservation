#! /usr/bin/perl
#
######################################################################################################################################################################
#### Program:     Reads sequences and removes duplicates, remove sequences with Ns and remove sequences with only 1nt. 
#### Input:       1. Fasta files (having a header and sequence in one-line. If not, run the script "fasta_oneLiner.sh" to convert multi-line fasta to one-line.
####Output:       1. Sequences in fasta format.
##
####Notes:        The program is executed as a job script on UBELIX server using the modules available on the server and hence, 
###               run the job batch script named:process_windows_job.sh to use this program.
####Author:       Mehta D. (dolly.mehta@unibe.ch)
####Date:         September 09, 2023
########################################################################################################################################################################

open(FA, $ARGV[0]);
@lines = <FA>;
close(FA);

my %windows={};
my @N_keys=();

for(my $i=0; $i<$#lines; $i+=2)
{
	$j=$i+1;
	chomp($lines[$i]);
	$windows{$lines[$i]}=$lines[$j];
}
undef(@lines);

#Delete any sequences with Ns
for my $key (keys %windows)
{
	push(@N_keys, $key) if($windows{$key} =~ m/^$/);
	push(@N_keys, $key) if($windows{$key} =~ m/n|N/);
}
delete @windows{@N_keys};

$intm1 = "$ARGV[1]"."_noNs.fasta";
open(N_FA, ">$intm1");

for my $key (keys %windows)
{
	print N_FA "$key\n$windows{$key}";
}

close(N_FA);

#Delete poly-A/T/G/C sequences
print "ID\tseqLen\tcount_A\tcount_T\tcount_G\tcount_C\tcount_U\n";

for my $key (keys %windows)
{
	my $count_A=$count_T=$count_G=$count_C=0;
	my $val = uc($windows{$key});
	my $len = length($val)-1;

	push(@polyKeys, $key) if($val =~ m/^$/);
	foreach my $s (split(//, $val))
	{
		$count_A+=1 if($s eq 'A');
		$count_T+=1 if($s eq 'T');
		$count_G+=1 if($s eq 'G');
		$count_C+=1 if($s eq 'C');
		$count_U+=1 if($s eq 'U');
	}

	if($count_A==$len || $count_T==$len || $count_G==$len || $count_C==$len || $count_U==$len)
	{
		push(@polyKeys, $key);
	}
	else
	{
		print "$key\t$len\t$count_A\t$count_T\t$count_G\t$count_C\t$count_U";
	}	
	print "\n";
}
delete @windows{@polyKeys};

$intm2 = "$ARGV[1]"."_noNs_polyN.fasta";
open(POLYN, ">$intm2");

for my $key (keys %windows)
{
	print POLYN "$key\n$windows{$key}";
}

close(POLYN);

undef(@polyKeys), undef(@N_keys);
