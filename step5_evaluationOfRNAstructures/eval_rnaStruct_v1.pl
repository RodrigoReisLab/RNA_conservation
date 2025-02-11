#! usr/bin/perl
use POSIX;

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
if($total_open_brackets>=5 && $total_open_brackets<10)
{
	if($ssConsensus =~ /($stem_open_brackets){$total_open_brackets}/)
	{
		#check for each candidate sequence in the cluster is having the base-pair
		my $threshold = $total_open_brackets;
		eval_struct(\@alignments, $threshold);
	}
	else
	{
		print "$cluster\t$nalns\t$ssConsensus\t$ss_consLen\t$total_open_brackets\t$passed_candidates_count\tLow\n";
	}
}
elsif($total_open_brackets>=10)
{
	#Note only checking for 50% of open brackets, since at SS_cons the total number of open and close brackets are identical
	my $result = $total_open_brackets * 0.5;
	my $threshold = $result >= int($result) + 0.5 ? ceil($result) : floor($result);
	eval_struct(\@alignments, $threshold);
}
else
{
	print "$cluster\t$nalns\t$ssConsensus\t$ss_consLen\t$total_open_brackets\t$passed_candidates_count\tLow\n";
}

#Write passed candidates to stockholm file
if($#passed_cands > 0)
{
	open(STO, ">$ARGV[1]");
	foreach my $t (@top_lines)
	{
	        $t =~ s/^\s//;
		$t= "#=GF SQ $passed_candidates_count\n" if($t =~ /SQ/);
	        print STO $t;
	}
	print STO "#=GF ID $cluster\n#=GF CC Processed Motif region (trimAlignment.pl) are screened and passed candidates are written \n\n";
	foreach $passed (@passed_cands)
	{
		print STO "$passed\n";
	}
	print STO "#=GC SS_cons\t$ssConsensus\n";
	print STO "//\n";

	close(STO);
}

# sub-routine to evaluate structure
sub eval_struct
{
	my ($aln_ref, $bp_threshold) = @_;
	@input_aln = @$aln_ref;

	#check candidate sequence having 50% base-pair
	open(OUT, ">$ARGV[2]");
	print OUT "seqID\tbp_threshold\tTotalBPs\t#OpenBPs\t#ClosingBPs\t#BPsActuallyForming\tEvaluation\n";
	$nucl = qr/[ATCGUatcgu]/;
	foreach $cand (@input_aln)
	{
		my ($seq_id, $seq_cand) = split(/\t/, $cand);
		my $open_brack_char = $close_brack_char = "";
		my @open_brack_chars = @close_brack_chars = @o_stem = @c_stem = ();
		for my $i (0 .. $#open_bracket_indices)
		{
			my $open_brack_loc = $open_bracket_indices[$i];
			my $close_brack_loc = $close_bracket_indices[$i];
			my $next_open = $open_bracket_indices[$i+1] if($i+1 <= $#open_bracket_indices);
			$next_open = $close_bracket_indices[$#open_bracket_indices] if($i+1 > $#open_bracket_indices);
			$open_brack_char = substr($seq_cand, $open_brack_loc, 1);
			$close_brack_char = substr($seq_cand, $close_brack_loc, 1);

			if($next_open < $close_brack_loc)
			{
				push(@o_stem, $open_brack_loc);
				push(@c_stem, $close_brack_loc);
			}
			else
			{
				#Getting the last base-pair of the stem
				push(@o_stem, $open_brack_loc);
				push(@c_stem, $close_brack_loc);
				#Checking the base-pair potential between the nucleotides, this will affect the base-pairs that each sequence can form.
				$effective_bp_count += eval_stem_loop(\@o_stem,\@c_stem, $seq_cand);
				undef(@o_stem), undef(@c_stem);

			}

			push(@o_nucl, $open_brack_char) if($open_brack_char =~ $nucl);
			push(@c_nucl, $close_brack_char) if($close_brack_char =~ $nucl);
			push(@open_brack_chars, $open_brack_char);
			push(@close_brack_chars, $close_brack_char);
		}

		my $o_nucl_count = scalar @o_nucl;
		my $c_nucl_count = scalar @c_nucl;

		my $nucl_count = $o_nucl_count < $c_nucl_count ? $o_nucl_count : $c_nucl_count;
		#Check for sequences having $threshold base-pairs and write to output.
		if($effective_bp_count>=$bp_threshold)
		{
			print OUT "$seq_id\t$bp_threshold\t$total_open_brackets\t$o_nucl_count\t$c_nucl_count\t$effective_bp_count\tPassed\n";
			#Write the sequence to stockholm file
			push(@passed_cands, $cand);
		}
		else
		{
			print OUT "$seq_id\t$bp_threshold\t$total_open_brackets\t$o_nucl_count\t$c_nucl_count\t$effective_bp_count\tFail\n";
		}
		undef(@open_brack_chars), undef(@close_brack_chars), undef(@o_nucl), undef(@c_nucl), undef(@motif), undef($effective_bp_count);
	}
	close(OUT);

	# Checking if the alignment has enough candidates to call it an alignment
	$passed_candidates_count = $#passed_cands + 1;
	print "$cluster\t$nalns\t$ssConsensus\t$ss_consLen\t$total_open_brackets\t$passed_candidates_count\tHigh\n" if($passed_candidates_count >= 10);
	print "$cluster\t$nalns\t$ssConsensus\t$ss_consLen\t$total_open_brackets\t$passed_candidates_count\tMid\n" if($passed_candidates_count>=7 && $passed_candidates_count < 10);
	print "$cluster\t$nalns\t$ssConsensus\t$ss_consLen\t$total_open_brackets\t$passed_candidates_count\tLow\n" if($passed_candidates_count < 7);

}

sub eval_stem_loop
{
	my ($open_ref, $close_ref, $seq) = @_;
	$bp_passed=0;
	for my $i (0 .. $#{$open_ref})
	{
		my $open_index = $open_ref->[$i];
		my $close_index = $close_ref->[$#{$close_ref} - $i];

		$open_char = substr($seq, $open_index, 1);
		$close_char = substr($seq, $close_index, 1);

		#Checking base-pair potential
		if(($open_char eq 'A' && $close_char eq 'U') || ($open_char eq 'U' && $close_char eq 'A') || ($open_char eq 'C' && $close_char eq 'G') || ($open_char eq 'G' && $close_char eq 'C') || ($open_char eq 'U' && $close_char eq 'G') || ($open_char eq 'G' && $close_char eq 'U'))
		{
			$bp_passed += 1;
		}
	}

	#Evaluating if each stem of the RNA has 75% structure
	my $result = $#{$open_ref} * 0.75;
	my $stem_thresh = $result >= int($result) + 0.5 ? ceil($result) : floor($result);

	$bp_passed = 0 if($bp_passed < $stem_thresh);
	if($bp_passed<3)
	{
		$bp_passed = 0;
	}
	elsif($bp_passed==3)
	{
		$open_last = $open_ref->[-1];
		$close_first = $close_ref->[0];
		my $loop_nucl = $close_first - ($open_last+1);
		my $loop=substr($seq, $open_last+1, $loop_nucl);
		$loop =~ s/-//g;
		$bp_passed = 0 if(length($loop) > 5);
	}
	return($bp_passed);

}
undef($ssConsensus), undef(@alignments), undef(@aln), undef(@open_bracket_indices), undef(@close_bracket_indices), undef(@passed_cands);





