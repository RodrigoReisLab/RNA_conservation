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

$input_file = $ARGV[0].".fasta";
open(FA, $input_file);

my %sequences;

my $header = "";
my $sequence = "";

while (my $line = <FA>)
{
	chomp $line;
	if ($line =~ /^>/)
	{
		# If a new header is encountered, process the previous sequence
		add_sequence();
		# Start processing a new header
		$header = $line;
	}
	else
	{
		# Concatenate sequence lines
		$sequence .= $line;
	}		
}
# Process the last sequence in the file
add_sequence();

close(FA);


my $out_file = $ARGV[0] . '_uniq.fasta';
open(OUT, ">$out_file") or die "Could not open $out_file: $!\n";

foreach my $seq (keys %sequences) {
    print OUT "$sequences{$seq}\n$seq\n";
}

undef(%sequences);
close(OUT);

sub add_sequence
{
	if ($header && $sequence)
	{
		if (!exists $sequences{$sequence}) 
		{
         		$sequences{$sequence} = $header;
        	}
		$sequence = "";  # Reset the sequence
	}	
}
