#! /usr/bin/perl

#####################################################################################################################################################################
### Program:     Take the genome sequence and create windows each of user-specified length with a overl 
### Input:       1. Fasta files (having a header and sequence in one-line. If not, run the script "fasta_oneLiner.sh" to convert multi-line fasta to one-line.
###Output:       1. Windowed sequences of length defined by user.
#
###Notes:        The program is executed as a job script on UBELIX server using the modules available on the server and hence, 
##               run the job batch script named: windowsCreater_job.sh to use this program.
###Author:       Mehta D. (dolly.mehta@unibe.ch)
###Date:         July 31, 2023
#######################################################################################################################################################################

use Getopt::Long;

my $fasta, $window, $overlap;

GetOptions( 
	'fasta|f=s'   => \(  $fasta = undef ),
	'window|w=s'    => \( $window = undef ),
	'overlap|p=s'    => \( $overlap = undef ),
	'outfile|O=s'	=> \( $out = undef ),
)
    or usage(); 

if (not defined $fasta or not defined $window or not defined $out or $#ARGV == 0)
{
	say STDERR "Input file, window length and output file are mandatory\n";
	usage();
}

my %windows={};
open(FA, $fasta);
@lines = <FA>;
close(FA);

open(OUT, ">$out");

for(my $i=0; $i<$#lines; $i+=2)
{
	$j=$i+1;
	chomp($lines[$i]);
	$windows{$lines[$i]}=$lines[$j];
}

for my $key (keys %windows)
{
	($k, $junk) = split(/ /, $key);
	my $counter=$slider=$repos=0;
	$seq=$windows{$key};
	$slider=$window;
	$size=length($seq);
	
	$overlap = 0 if(not defined $overlap); #Windows with no overlap
	for(my $start=0; ($repos+$window)<$size; $start=$repos)
	{
		my $start_coord=$repos;
		my $end_coord=$start_coord + $window;
		my $header = ${k}."_".$counter."_".$start_coord."_".$end_coord;
		if($slider<$size)
		{
			$slider=$start+$window;
			$windowed_seq = substr( $seq, $repos, $window);
			$repos=$slider-$overlap;
			print OUT "$header\n$windowed_seq\n";
			undef($windowed_seq);
		}
		$counter+=1;
	}
}

close(OUT);
undef(@lines), undef(%windows), undef($overlap); 

sub usage {
	say STDERR "Usage: $0 -f inputFastaFile -w window_length -p overlap_length [optional] -O output_filename \nPlease maintain the same order of inputs\n";   # full usage message
	exit;
}
