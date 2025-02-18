#! usr/bin/perl

#Reads 1 input alignment file (stockholm format) and outputs 3 files (cleaned alignment: stockholm, evaluation of each RNA candidate and the total cleaned alignment: tab-separated)

#Getting cluster name
my $cluster = $ARGV[0];
$cluster =~ s/_motif\.sto$//;;

open(ALN, $ARGV[0]);
@aln = <ALN>;
close(ALN);

#Getting top header lines from the input file. Mandatory for maintaining the format for output file.
$n=2; #Modify according to your requirement. Since PERL is 0 based, we take total number -1.
my @top_lines = @aln[0..$n];

#Get secondary structure of the RNA family (ss_cons line in stockholm files)
my @ss_cons = map { "$_\t$aln[$_]" } grep { $aln[$_] =~/\bSS_cons\b/i } 1 .. @aln;
$ss_cons_len = @ss_cons;

#Check if the alignment is wrapped. It will be typically if the output is from mlocarna. This is important to check to calculate the motif
my $motif_loc="";
my $ssConsensus=$ss_consLen="";
my $stem_open_brackets = qr/[\(\{\[<]/;
my $stem_close_brackets = qr/[\)\}\]>]/;

if($ss_cons_len == 1)
{
	$ss_cons[0] =~ s/\t{2,}|\s+/\t/g;
	my ($lineIdx, $word, $w1, $p) = split(/\t/, $ss_cons[0]);
	$ssConsensus = $p;
	$ss_consLen = length($ssConsensus);

	#Get the motif start position in the alignment. By default we check for first open bracket: ( or { or [ or < - as an indication of a stem.
	my $stem_start_pos = index($ssConsensus, $1) if $ssConsensus =~ /($stem_open_brackets)/;
	my $stem_end_pos = rindex($ssConsensus, $1) if $ssConsensus =~ /($stem_close_brackets)/;

	$motif_loc = "$stem_start_pos\t$stem_end_pos";
}
else
{
	#Concatenate wrapped alignments in a line
	@sorted = sort { $a <=> $b } @ss_cons; # Extra precautionary step to making sure lines are in the order, else it will goof up the alignment
	foreach $line (@sorted)
	{
		$line =~ s/\t{2,}|\s+/\t/g;
		my ($lineIdx, $word, $w1, $p) = split(/\t/, $line);
		$ssConsensus = $ssConsensus.$p;
	}
	$ss_consLen = length($ssConsensus); 

	#Get the motif start position in the alignment. By default we check for first open bracket: ( or { or [ or < - as an indication of a stem.
	 my $stem_start_pos = index($ssConsensus, $1) if $ssConsensus =~ /($stem_open_brackets)/;
	 my $stem_end_pos = rindex($ssConsensus, $1) if $ssConsensus =~ /($stem_close_brackets)/;

	 $motif_loc = "$stem_start_pos\t$stem_end_pos";
}
chomp($ssConsensus);

#Get the IDs of the sequences in the alignment
foreach $line (@aln)
{
	if($line !~ /^#|^\//)
	{
		$count +=1;
		next if($line !~ /[ATUGCatugc]/);
		$line =~ s/\t{2,}|\s+/\t/g;
		my ($id,$seq) = split(/\t/, $line) if($line ne "");
		$entry = $count."\t".$id."\t".$seq;
		push(@aln_entries, $entry) if($line =~ /[ATUGCatugc]/);
		
		#For wrapped alignment
		push(@ids, $id) if(!grep( /\b$id\b/, @ids));
	}
}
my $nalns = @aln_entries;
my $nids = @ids;

#Get the sequence concatenated for wrapped alignments, retaining the alignment, while for single-line alignment - get the sequences as aligned (Meaning, gaps are retained).
my @alignments=();
if($nids != $nalns)
{
	foreach my $key (@ids)
	{
		my $concat_seq = "";
		my @broken_seqs = grep { /\b$key\b/ } @aln_entries;
		$srNo=0, $id="", $p="";
		foreach $piece (@sorted_brokenSeqs)
		{
			($srNo, $id, $p) = split(/\t/, $piece);
			$concat_seq = $concat_seq.$p;
		}
		my $entry = $id."\t".$concat_seq;
		push(@alignments, $entry);
	}
}
else
{
	foreach my $candidates (@aln_entries)
	{
		my ($srNo, $id, $p) = split(/\t/, $candidates);
		my $entry = $id."\t".$p;
		push(@alignments, $entry);
	}
}

#Check the alignments for following characteristics:
#1. The motif structure contains a total minimum 5 base-pairs. If the bp_count is 5, than the base-pairs should be consective.
#2. Each sequence in the alignment is evaluated for the structure. if it is 5bp, sequences should contain the entire motif,  else the sequences should contain 50% of the structure
#3. The number of sequences at the end of evaluation should be atleast 7 sequences. More the sequences, merrier it is!

#Counting number of opening base-pairs
my @open_bracket_indices, @close_bracket_indices;
while ($ssConsensus =~ /$stem_open_brackets/g)
{
	my $index = $-[0];
	push @open_bracket_indices, $index;
}
while ($ssConsensus =~ /$stem_close_brackets/g)
{
	my $index = $-[0];
	push @close_bracket_indices, $index;
}
my $total_open_brackets = @open_bracket_indices;
my $total_close_brackets = @close_bracket_indices;
my @passed_cand;
my $passed_candidates_count = 0;

#Evaluating structure
print "clusterFile\tNseqs\tssConsensus\tss_consLen\tTotal_basepairs\tNseqsPassingEval\tConfidenceOnStructure\n";
if($total_open_brackets==5)
{
	if($ssConsensus =~ /($stem_open_brackets){5}/)
	{
		#check for each candidate sequence in the cluster is having the base-pair
		my $threshold = $total_open_brackets;
		eval_struct(\@alignments, $threshold);
	}
}
elsif($total_open_brackets>5)
{
	#Note only checking for 50% of open brackets, since at SS_cons the total number of open and close brackets are identical
	my $threshold = $total_open_brackets * 0.5;
	eval_struct(\@alignments, $threshold);
}
else
{	
	$passed_candidates_count = 0;
	print "$cluster\t$nalns\t$ssConsensus\t$ss_consLen\t$total_open_brackets\t$passed_candidates_count\tLow\n";
}

#Write passed candidates to stockholm file
open(STO, ">$ARGV[1]");
foreach my $t (@top_lines)
{
        $t =~ s/^\s//;
	$t= "#=GF SQ $passed_candidates_count\n" if($t =~ /SQ/);
        print STO $t;
}
print STO "#=GF Processed Motif region (trimAlignment.pl) are screened and passed candidates are written \n\n";
foreach $passed (@passed_cands)
{
	print STO "$passed\n";
}
print STO "#=GC SS_cons\t$ssConsensus\n";
print STO "//\n";

close(STO);

# sub-routine to evaluate structure
sub eval_struct
{
	my ($aln_ref, $bp_threshold) = @_;
	@input_aln = @$aln_ref;

	#check candidate sequence having 50% base-pair
	open(OUT, ">$ARGV[2]");
	print OUT "seqID\tbp_threshold\tTotalBPs\t#OpenBPs\t#ClosingBPs\tEvaluation\n";
	foreach $cand (@input_aln)
	{
		my ($seq_id, $seq_cand) = split(/\t/, $cand);
		my $open_brack_char=my $close_brack_char="";
		my @open_brack_chars, @close_brack_chars;
		for my $i (0 .. $#open_bracket_indices)
		{
			my $open_brack_loc = $open_bracket_indices[$i];
			my $close_brack_loc = $close_bracket_indices[$i];
			$open_brack_char = substr($seq_cand, $open_brack_loc, 1);
			$close_brack_char = substr($seq_cand, $close_brack_loc, 1);

			push(@open_brack_chars, $open_brack_char);
			push(@close_brack_chars, $close_brack_char);
		}
		$nucl = qr/[ATCGUatcgu]/;
		my @o_nucl = grep { $_ =~ $nucl } @open_brack_chars;
		my @c_nucl = grep { $_ =~ $nucl } @close_brack_chars;

		my $o_nucl_count = scalar @o_nucl;
		my $c_nucl_count = scalar @c_nucl;

		#Check for sequences having $threshold base-pairs and write to output.
		if($o_nucl_count>=$bp_threshold and $c_nucl_count>=$bp_threshold)
		{
			print OUT "$seq_id\t$bp_threshold\t$total_open_brackets\t$o_nucl_count\t$c_nucl_count\tPassed\n";
			#Write the sequence to stockholm file
			push(@passed_cands, $cand);
		}
		else
		{
			print OUT "$seq_id\t$bp_threshold\t$total_open_brackets\t$o_nucl_count\t$c_nucl_count\tFail\n";
		}
		undef(@open_brack_chars), undef(@close_brack_chars), undef(@o_nucl), undef(@c_nucl);
	}
	close(OUT);

	# Checking if the alignment has enough candidates to call it an alignment
	$passed_candidates_count = $#passed_cands + 1;
	print "$cluster\t$nalns\t$ssConsensus\t$ss_consLen\t$total_open_brackets\t$passed_candidates_count\tHigh\n" if($passed_candidates_count >= 10);
	print "$cluster\t$nalns\t$ssConsensus\t$ss_consLen\t$total_open_brackets\t$passed_candidates_count\tMid\n" if($passed_candidates_count < 10);

}
undef($ssConsensus), undef(@alignments), undef(@aln), undef(@open_bracket_indices), undef(@close_bracket_indices), undef(@passed_cands);





