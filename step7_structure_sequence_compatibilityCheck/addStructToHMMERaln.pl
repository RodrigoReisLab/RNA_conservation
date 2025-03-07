#! /usr/bin/perl

#one sequence stockholm (obtained from esl-alimask command). These entries must not be wrapped. It should be in one line.
open(ONE, $ARGV[0]);
@one = <ONE>;
close(ONE);

@one_cons_arr = grep { /#=GC SS_cons/ } @one;
my $temp = $one_cons_arr[0];
$temp =~ s/\s\s\s*/\t/g;
my ($m, $one_cons) = split(/\t/, $temp);

#alignment of RNAs obtained after nhmmer (hmmsearch) run in stockholm format (file usually named file_nhmmer.sto'). Can be multi-line as well!
open(STO, $ARGV[1]);
@sto = <STO>;
@rf_line = grep { /#=GC RF/ } @sto;

my $prev_xcount = my $line_counter = 0;
my %cons;
for my $rf (@rf_line) {
	my $sub_cons;
	$rf =~ s/\s\s\s*/\t/g;
	my ($id, $line) = split(/\t/, $rf);
	chomp($line);

	my $x_count = () = $line =~ /x/g;

	#Writing consensus guided by RF line
	my @rf_chars = split(//, $line);
	my $counter=$prev_xcount;
	for my $char (@rf_chars)
	{
		my $cons_char;
		if($char eq 'x') {
			$cons_char = substr($one_cons, $counter, 1);
			$counter++;
		}
		else
		{
			$cons_char='.';
		}
		$sub_cons .= $cons_char;
	}
	$prev_xcount = $counter;
	#print "$line_counter\t$sub_cons\n";
	$cons{$line_counter} = $sub_cons;
	$line_counter++;

	undef($sub_cons), undef(@rf_chars);
}
undef($prev_xcount);

#Adding cons line to the sto
my $cons_counter = 0;
my @mod_sto;
foreach my $sto_line (@sto)
{
	push(@mod_sto, $sto_line);
	if($sto_line =~ /#=GC RF/) {
		if(exists $cons{$cons_counter}) {
			$cons_line = "#=GC SS_cons                               $cons{$cons_counter}";
			push(@mod_sto, $cons_line);
			$cons_counter++;
		}
	}
}

foreach my $l (@mod_sto) {
	chomp($l);
	print "$l\n";
}

undef(%cons), undef(@sto), undef(@rf), undef(@one), undef(@mod_sto), undef(@rf_line), undef(@one_cons_arr);
close(STO);
