#! usr/bin/perl

#input list of clusters 
open(FA, $ARGV[0]);
my @fa = <FA>;
close(FA);

#write output file
open(OUT, ">$ARGV[1]");

my $len = @fa;

# Removing gaps from the FASTA (Hint: Input files are obtained as input from trimAlignment.pl)
for(my $i=0; $i < $len; $i=$i+2)
{
	my $header=$fa[$i];
	chomp($header);
	my $a=$i+1;
	my $seq=$fa[$a];
	$seq =~ s/-//g;
	chomp($seq);
	print OUT "$header\n$seq\n";
}

close(OUT);
