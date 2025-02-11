#! usr/bin/perl
use POSIX;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

#Reads 1 input alignment file (stockholm format) and outputs 3 files (cleaned alignment: stockholm, evaluation of each RNA candidate and the total cleaned alignment: tab-separated)

my ($ext, $motif_thresh, $bp_thresh, $gc_bp_thresh, $help, $man);
#Getting command line options
GetOptions(
	'e|extension=s'		=> \$ext, 		#any flanking in filename that should be replaced. mandatory to remove .sto from extension - to obtain clean output filename
	'mt|motif_threshold=f'	=> \$motif_thresh, 	#total motif threshold [0.0-1.0] => set to 0.5 for evaluation/cleaning step
	't|bp_threshold=f'	=> \$bp_thresh, 	#each stem in the motif should have x% of base-pairs [0.0-1.0] => set to 0.75 for evaluation/cleaning step
	'gc|bp_gc=f'		=> \$gc_bp_thresh, 	#each stem in the motif should have x% of GC/CG base-pairs [0.0-1.0] -> set to 0.30
	'd|duplication=i'	=> \$dupl_flag,		#A flag to include duplications in the alignment and all the calculations: 0: duplications OFF - removes duplicationsand 1: duplications ON, retains duplications in the alignment
	'h|help'		=> \$help,
	'manual'		=> \$man
 ) or die "Usage: Incorrect Usage.\n Use --help or --manual to see the correct use of the script\n";

#Display help if needed
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

# Check if the mandatory extension option is provided
if(!defined $ARGV[0])
{
	die "Error: Input alignment file not given\n";
}
if (!defined $ext)
{
	die "Error: The extension option for filename is mandatory.\nExample: -f is 'a_motif.sto'. Set -e to '_motif.sto'\n";
}
if (!defined $dupl_flag)
{
	die "Error: Setting duplication flag is mandatory.\nSet -d to 0 (remove duplicate sequences) or 1 (retain duplications)\n";
}
if (!defined $motif_thresh)
{
	print "Motif threshold not provided. Setting motif_threshold to default of 0.5.\n";
	$motif_thresh = 0.5;
}
if (!defined $bp_thresh)
{
	print "Base-pairs threshold for stems in motif is not provided. Setting stem_bp_threshold to default of 0.75.\n";
	$bp_thresh = 0.75;
}
if (!defined $gc_bp_thresh)
{
		print "%GC/CG threshold in base-pairs is not provided. Setting gc_bp_threshold to default of 0.30.\n";
		$gc_bp_thresh = 0.30;
}

open(ALN, $ARGV[0]);
@aln = <ALN>;
close(ALN);

#Getting cluster name
my $cluster = $ARGV[0];
$cluster =~ s/${ext}$//;
$cluster =~ s/^allchloroGenomes_// if($cluster =~ m/^allchloroGenomes_/);

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
my @alignments;
if($nids != $nalns)
{
	foreach my $key (@ids)
	{
		my $concat_seq = "";
		my @broken_seqs = grep { /\b$key\b/ } @aln_entries;
		my @sorted_brokenSeqs = sort { $a <=> $b } @broken_seqs;
		my $srNo=0, $id="", $p="";
		foreach $piece (@sorted_brokenSeqs)
		{
			($srNo, $id, $p) = split(/\t/, $piece);
			$concat_seq = $concat_seq.$p;
		}
		my $entry = $id."\t".$concat_seq;
		push(@alignments, $entry);
		undef(@broken_seqs), undef(@sorted_brokenSeqs);
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

#Removing identical motif sequences (retaining one as a representative - typically first that is seen in the alignment)
my @uniq_aln, @duplRNAcands;
if($dupl_flag == 0)
{
	my ($uniq_ref, $dupl_ref) = motif_seen(\@alignments);
	@uniq_aln = @$uniq_ref;
	@duplRNAcands = @$dupl_ref
}
else
{
	@uniq_aln = @alignments;
}
my $uniqAln_size = @uniq_aln;

#Check the alignments for following characteristics:
#1. The motif structure contains a total minimum 5 base-pairs. If the bp_count is 5, than the base-pairs should be consective.
#2. Each sequence in the alignment is evaluated for the structure. if it is 5-9bp, sequences should contain the entire motif,  else the sequences should contain 50% of the structure
#3. The number of sequences at the end of evaluation should be atleast 7 sequences. More the sequences, merrier it is!

#Extracting hairpins step-by-step
#Extracting base-pairs
my $pairs = extract_pairs_from_cons($ssConsensus);

#Making pair matrix and finding hairpins from the matrix. These may include hairpins which are nested and declared as 1 hairpin. Need to separate the hairpins and hence nesting hairpins need to be identified.
my $base_pair_matrix = create_pairMatrix($pairs);
open(my $mat, ">${cluster}_pair_matrix.txt") or die "Could not open file '${cluster}_pair_matrix.txt' $!";
foreach my $row (@$base_pair_matrix) {
	my $line = join("\t", @$row);  # Tab-separated values
	print $mat "$line\n";
}
close($mat);

# Find the outermost hairpin/ nested hairpins
my $final_hairpins = find_hairpin_with_pair($pairs);
open(my $hp, ">${cluster}_detected_hairpins.txt") or die "Could not open file '${cluster}_detected_hairpins.txt' $!";
foreach my $hairpin (@$final_hairpins) {
	my $line = join("\t", map { join(",", @$_) } @$hairpin);
	print $hp "$line\n";
}
close($hp);

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
print "clusterFile\tNseqs\tNuniqSeqs\tssConsensus\tss_consLen\tTotal_basepairs\tNseqsPassingEval\tConfidenceOnStructure\n";
if($total_open_brackets>=5 && $total_open_brackets<10)
{
	if($ssConsensus =~ /($stem_open_brackets){$total_open_brackets}/)
	{
		#check for each candidate sequence in the cluster is having the base-pair
		my $threshold = $total_open_brackets;
		eval_struct(\@uniq_aln, $threshold);
	}
	else
	{
		print "$cluster\t$nids\t$uniqAln_size\t$ssConsensus\t$ss_consLen\t$total_open_brackets\t$passed_candidates_count\tLow - No evaluation performed as number of sequences<5 for evaluation\n";
	}
}
elsif($total_open_brackets>=10)
{
	#Note only checking for 50% of open brackets, since at SS_cons the number of open and close brackets are identical
	my $result = $total_open_brackets * $motif_thresh;
	my $threshold = $result >= int($result) + 0.5 ? ceil($result) : floor($result);
	eval_struct(\@uniq_aln, $threshold);
}
else
{
	print "$cluster\t$nids\t$uniqAln_size\t$ssConsensus\t$ss_consLen\t$total_open_brackets\t$passed_candidates_count\tLow\n";
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

# sub-routines
# Removing duplicates
sub motif_seen
{
	my $input_aln_ref = $_[0];
	my @input_alns = @$input_aln_ref;
	my %motif_cands, %uniq_cands, %duplHits, @duplRNAcands;

	foreach $in (@input_alns)
	{
		my ($key_seqID, $value_seq) = split(/\t/, $in);
		$motif_cands{$key_seqID} = $value_seq;
	}

	foreach my $k (keys %motif_cands)
	{
		my $motif_seq = $motif_cands{$k};

		#Removing duplicate motif sequences
		if (!grep { $_ eq $motif_seq } values %uniq_cands)
		{
			$uniq_cands{$k} = $motif_seq;
		}
		else
		{
			#Getting coordinates for RNAs with identical sequence homologs
			foreach my $uniq_key (keys %uniq_cands)
			{
	push(@{ $duplHits{$uniq_key} }, $k) if($uniq_cands{$uniq_key} eq $motif_seq)
			}
		}

	}
	
	foreach my $key (keys %uniq_cands)
	{
		my $uniq_cand_entry = $key."\t".$uniq_cands{$key};
		push(@uniq_arr, $uniq_cand_entry);
	}

	foreach my $key (keys %duplHits)
	{
		my $cands = join(",", @{ $duplHits{$key} });
		my $identRNAcands=$key."\t".$cands;
		push(@duplRNAcands, $identRNAcands);
	}

	undef(@input_alns), undef(%motif_cands), undef(%uniq_cands), undef(%duplHits);
	return(\@uniq_arr, \@duplRNAcands);
}

#Getting base-pairs and identify hairpins from consensus structure
sub extract_pairs_from_cons {
		my $structure = $_[0];
		my @stack, @pairs, @tmp;
		for my $i (0 .. length($structure) - 1)
		{
	my $char = substr($structure, $i, 1);
	if ($char =~ /$stem_open_brackets/) 
	{
		push(@stack, $i);
	} 
	elsif ($char =~ /$stem_close_brackets/) 
	{
		my $hp_start = pop @stack;
		push(@tmp, [$hp_start, $i]);
	}
		}
		
		@pairs = sort { $a->[0] <=> $b->[0] } @{tmp};
	return \@pairs;
}

#Subroutine to find the hairpin containing the nesting boundary base-pair
sub find_hairpin_with_pair {
	my ($pairs_ref) = @_;
	my @pairs = @$pairs_ref;
	my @separated_hairpins;
	my %visited;

	# Construct the adjacency list for all pairs
	my %graph;
	foreach my $pair (@pairs) {
		my ($row, $col) = @$pair;
		$graph{"$row,$col"} = [];
	}

	# Build the graph (connect pairs that are adjacent based on are_connected)
	foreach my $pair (@pairs) {
		my ($row1, $col1) = @$pair;
		foreach my $neighbor (@pairs) {
			my ($row2, $col2) = @$neighbor;
			next if $pair eq $neighbor;
			if (are_connected($row1, $col1, $row2, $col2)) {
	push @{$graph{"$row1,$col1"}}, "$row2,$col2";
			}
		}
	}

	# Perform DFS to find all connected components
	sub dfs {
		my ($node, $component) = @_;
		return if $visited{$node};
		$visited{$node} = 1;
		push @$component, $node;
		foreach my $neighbor (@{$graph{$node}}) {
			dfs($neighbor, $component);
		}
	}

	my @components;
	foreach my $node (keys %graph) {
		next if $visited{$node};
		my @component;
		dfs($node, \@component);
		push @components, \@component;
	}

	# Convert components back to pairs and add to separated_hairpins
	foreach my $component (@components) {
		my @hairpin = map { [ split /,/, $_ ] } @$component;
		push @separated_hairpins, \@hairpin;
	}

	# Sort each hairpin
	foreach my $hairpin (@separated_hairpins) {
		@$hairpin = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @$hairpin;
	}

	# Sort all hairpins
	@separated_hairpins = sort {
		my $cmp = join(',', map { sprintf("%04d%04d", $_->[0], $_->[1]) } @$a);
		my $cmp2 = join(',', map { sprintf("%04d%04d", $_->[0], $_->[1]) } @$b);
		$cmp cmp $cmp2;
	} @separated_hairpins;
	
	#print ">>>> Separated hairpins: \n", Dumper(\@separated_hairpins);
	return \@separated_hairpins;
}

# Determine if two pairs are connected (more complex criteria)
sub are_connected {
	my ($row1, $col1, $row2, $col2) = @_;
	# Flexible connectivity criteria, based on relative positions
	return (abs($row1 - $row2) <= 12 && abs($col1 - $col2) <= 12) ||  # Adjacent
		(abs($row1 - $row2) <= 12 && abs($col1 - $col2) <= 12) ||  # Diagonal
		(abs($row1 - $row2) <= 12 && abs($col1 - $col2) <= 12);	# L-shaped
}

sub create_pairMatrix {
	my ($pairs) = @_;
	
	# Determine the matrix size
	my $max = 0;
	foreach my $pair (@$pairs) {
		$max = $pair->[1] if $pair->[1] > $max;
	}
	$max += 10;

	# Initialize the matrix
	my @matrix;
	for my $i (0..$max) {
		for my $j (1..$max) {
			$matrix[$i][$j] = 0;
		}
	}

	# Fill the matrix with the provided pairs
	foreach my $pair (@$pairs) {
		$matrix[$pair->[0]][$pair->[1]] = 1;
	}
	
	return \@matrix;
}

#Evaluating structures
sub eval_struct
{
	my ($aln_ref, $bp_threshold) = @_;
	@input_aln = @$aln_ref;

	#check candidate sequence having 50% base-pair of the total structure
	open(OUT, ">$ARGV[2]");
	print OUT "seqID\tlength\tGCcontent\tbp_threshold\tTotalBPs\t#OpenBPs\t#ClosingBPs\tHairpinBPs_per\t#BPsActuallyForming\tRedundantRNAcandidates\tEvaluation\n";
	$nucl = qr/[ATCGUatcgu]/;
	
	#Evaluating each candidate in the alignment (removed/retained identical replicates)
	foreach $cand (@input_aln)
	{
		my ($seq_id, $seq_cand) = split(/\t/, $cand);
		my $seq_gc = gcContent($seq_cand);
		my $temp = $seq_cand;
		$temp =~ s/[\.,:_;-]//g;
		my $seq_len = length($temp);
		my $open_brack_char = $close_brack_char = "";
		my $total_bp_percent, $effective_bp_count;
		my @open_brack_chars = @close_brack_chars = @o_stem = @c_stem = ();
		for my $i (0 .. $#open_bracket_indices)
		{
			my $open_brack_loc = $open_bracket_indices[$i];
			my $close_brack_loc = $close_bracket_indices[$i];
			my $next_open = $open_bracket_indices[$i+1] if($i+1 <= $#open_bracket_indices);
			$next_open = $close_bracket_indices[$#open_bracket_indices] if($i+1 > $#open_bracket_indices);
			$open_brack_char = substr($seq_cand, $open_brack_loc, 1);
			$close_brack_char = substr($seq_cand, $close_brack_loc, 1);
	
			push(@o_nucl, $open_brack_char) if($open_brack_char =~ $nucl);
			push(@c_nucl, $close_brack_char) if($close_brack_char =~ $nucl);
			push(@open_brack_chars, $open_brack_char);
			push(@close_brack_chars, $close_brack_char);

		}
		#Getting each stem-loop (hairpin) of the motif for each RNA candidates
		foreach my $hp (@$final_hairpins)
		{
			foreach my $pair (@$hp) 
			{
				my ($open, $close) = @$pair;
				push(@o_stem, $open);
				push(@c_stem, $close);
			}
			#Checking the base-pair potential between the nucleotides, this will affect the base-pairs that each sequence can form.
			my ($sl_bp_count, $sl_bp_percent) = eval_stem_loop($hp, $seq_cand);
			$effective_bp_count += $sl_bp_count;
			$total_bp_percent .= $sl_bp_percent;
			$total_bp_percent .= "," unless $i == $#open_bracket_indices;
			undef(@o_stem), undef(@c_stem), undef($sl_bp_count);
		}

		my $o_nucl_count = scalar @o_nucl;
		my $c_nucl_count = scalar @c_nucl;

		my $nucl_count = $o_nucl_count < $c_nucl_count ? $o_nucl_count : $c_nucl_count;

		#Getting list of redundant RNA candidates if present
		my $redundantRNAcands="";
		my @parts;
		my ($duplicates) = grep { /$seq_id/ } @duplRNAcands;
		my @parts = split(/\t/, $duplicates);
		shift @parts;
		$redundantRNAcands = join(',', @parts);
		undef(@parts), undef($part);

		#Check for sequences having $threshold base-pairs and write to output.
		if($effective_bp_count>=$bp_threshold && $seq_gc > 0)
		{
			print OUT "$seq_id\t$seq_len\t$seq_gc\t$bp_threshold\t$total_open_brackets\t$o_nucl_count\t$c_nucl_count\t$total_bp_percent\t$effective_bp_count\t$redundantRNAcands\tPassed\n";
			#Write the sequence to stockholm file
			push(@passed_cands, $cand);
		}
		else
		{
			print OUT "$seq_id\t$seq_len\t$seq_gc\t$bp_threshold\t$total_open_brackets\t$o_nucl_count\t$c_nucl_count\t$total_bp_percent\t$effective_bp_count\t$redundantRNAcands\tFail\n";
		}
		undef(@open_brack_chars), undef(@close_brack_chars), undef(@o_nucl), undef(@c_nucl), undef(@motif), undef($effective_bp_count), undef($total_bp_percent);
	}
	close(OUT);

	# Checking if the alignment has enough candidates to call it an alignment
	$passed_candidates_count = $#passed_cands + 1;
	print "$cluster\t$nids\t$uniqAln_size\t$ssConsensus\t$ss_consLen\t$total_open_brackets\t$passed_candidates_count\tHigh\n" if($passed_candidates_count >= 10);
	print "$cluster\t$nids\t$uniqAln_size\t$ssConsensus\t$ss_consLen\t$total_open_brackets\t$passed_candidates_count\tMid\n" if($passed_candidates_count>=7 && $passed_candidates_count < 10);
	print "$cluster\t$nids\t$uniqAln_size\t$ssConsensus\t$ss_consLen\t$total_open_brackets\t$passed_candidates_count\tLow\n" if($passed_candidates_count < 7);

}

sub gcContent
{
	my $seq = $_[0];
	my $gc_count=0;
	
	while($seq =~ /G|C/ig)
	{
		$gc_count += 1;
	}
	
	return($gc_count);
}

sub eval_stem_loop
{
	#Check if it is a hairpin and no nested base-pairing is present
	my ($hp_ref, $seq) = @_;
	my $bp_passed=$gc=0;

	my @hp = @$hp_ref;
	#print Dumper \@hp;
	my $number_of_pairs = $#hp + 1;
	foreach my $pair (@hp)
	{
		my ($open_index, $close_index) = @$pair;
		$open_char = substr($seq, $open_index, 1);
		$close_char = substr($seq, $close_index, 1);

		#Checking base-pair potential
		if(($open_char eq 'A' && $close_char eq 'U') || ($open_char eq 'U' && $close_char eq 'A') || ($open_char eq 'C' && $close_char eq 'G') || ($open_char eq 'G' && $close_char eq 'C') || ($open_char eq 'U' && $close_char eq 'G') || ($open_char eq 'G' && $close_char eq 'U'))
		{
			$bp_passed += 1;

			if(($open_char eq 'C' && $close_char eq 'G') || ($open_char eq 'G' && $close_char eq 'C'))
			{
	$gc += 1;
			}
		}
	}
	#print "cand: $number_of_pairs\t$bp_passed\t$gc\n";
	#Evaluating if each stem of the RNA has 75% structure
	#Getting the minimum number of base-pairs needed for it to be called a stem.
	my $result = $number_of_pairs * $bp_thresh;
	my $stem_thresh = $result >= int($result) + 0.5 ? ceil($result) : floor($result);
	
	#Total number of base-pairs
	my $open_ref_length = $number_of_pairs;
	my $bp_passed_per = ($bp_passed / $open_ref_length) * 100;
	my $bp_passed_percent = sprintf("%.2f", $bp_passed_per);

	my $last_pair = $pairs[-1];
	$open_last = $last_pair->[0];
	$close_first = $last_pair->[1];
	my $loop_nucl = $close_first - ($open_last+1);
	my $loop=substr($seq, $open_last+1, $loop_nucl);
	$loop =~ s/[-\._~]//g;
	$bp_passed = 0 if($bp_passed < $stem_thresh);
	if($bp_passed<3)
	{
		$bp_passed = 0;
	}
	if($bp_passed==3)
	{
		$bp_passed = 0 if(length($loop) > 8);
	}

	#For all stems loop has to be atleast 2nt.
	$bp_passed = 0 if(length($loop) < 2);

	#Evaluating if each stem has 30% GC/CG base-pair, else the stems are weak
	my $bp_passed_for_gc = $bp_passed * $gc_bp_thresh;
	my $gc_thresh = $bp_passed_for_gc >= int($bp_passed_for_gc) + 0.5 ? ceil($bp_passed_for_gc) : floor($bp_passed_for_gc);

	#If GC base-pairs not part of stem, stems do not contribute to the motf, thus setting base-pairs as 0
	$bp_passed = 0 if($gc < $gc_thresh);
	
	return($bp_passed, $bp_passed_percent);

}
undef($ext), undef($motif_thresh), undef($bp_thresh), undef($gc_bp_thresh), undef($help), undef($man), undef($dupl_flag);
undef(@aln), undef($cluster), undef($n), undef(@top_lines), undef(@ss_cons), undef($ss_cons_len);
undef($motif_loc), undef($ssConsensus), undef($ss_consLen), undef($stem_open_brackets), undef($stem_close_brackets);
undef(@sorted), undef($line), undef($lineIdx), undef($word), undef($w1), undef($p), undef($stem_start_pos), undef($stem_end_pos);
undef(@ids), undef($count), undef($entry), undef(@aln_entries), undef($nalns), undef($nids), undef(@alignments), undef(@uniq_aln);
undef(@duplRNAcands), undef($uniqAln_size), undef(@passed_cand), undef($passed_candidates_count);
undef(@open_bracket_indices), undef(@close_bracket_indices), undef($total_open_brackets), undef($total_close_brackets);
undef(@uniq_arr), undef($final_hairpins), undef($concat_seq), undef(@broken_seqs);
undef(@sorted_brokenSeqs), undef($srNo), undef($id), undef($motif_seq), undef(%motif_cands), undef(%uniq_cands), undef(%duplHits);
undef(@uniq_arr), undef(@duplRNAcands), undef(@sorted_pairs), undef($nesting_pair), undef(@nesting_pairs), undef(%seen);
undef($structure), undef(@stack), undef(@pairs), undef(@tmp), undef(@components), undef(%visited), undef(%graph);

__END__

=head1 NAME

perl eval_rnaStruct_v2.pl - Reads an input alignment file in stockholm format, compares each sequence/ hit in the alignment with the consensus structure and evaluates if the sequence/hit can form the secondary structure. Evaluates at each hairpin (stem-loop), the total motif and the entire alignment. The estimates used for evaluation are output to evaluation.tsv (tab-separated) file for each cluster, a cleaned alignment (removing redundant and failed sequeces/ hits) and overall cluster confidence of alignment and other features is written to standard output.

=head1 SYNOPSIS

perl eval_rnaStruct_v2.pl [options] infile.sto outfile_cleaned.sto outfile_evaluation.tsv

 Options:
   -e, --extension		 Mandatory: any flanking in filename that should be replaced. (e.g., remove .sto)
   --mt, --motif_threshold   Total motif threshold [0.0-1.0]. Default: 0.5
   -t, --bp_threshold	  Each stem in the motif should have x% of base-pairs [0.0-1.0]. Default: 0.75
   --gc, --bp_gc		Each stem in the motif should have x% of GC/CG base-pairs [0.0-1.0]. Default: 0.30
   -d, --duplication	A flag to include duplications in the alignment and all the calculations: 0: duplications OFF - removes duplications and 1: duplications ON, retains duplications in the alignment
   --help		Display help message
   --manual	Display full manual

=head1 DESCRIPTION

B<This program> will read the given input file(s) and process them as described before. Expect 3 outputs:\n
1. evaluation.tsv	=> tab-seperated file with a header.
2. cleaned_aln.sto	=> alignment in stockholm format.
3. STDIN		=> cluster confidence for structure based on the alignment.

