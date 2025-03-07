#! usr/bin/perl

open(ALN, $ARGV[0]);
@aln = <ALN>;
close(ALN);

#Get header from the stockholm file (In my files it is the first 4 lines).
$n=3; #Modify according to your requirement. Since PERL is 0 based, we take total number -1.
my @top_lines = @aln[0..$n];

#Get secondary structure of the RNA family (ss_cons line in stockholm files)
my @ss_cons = map { "$_\t$aln[$_]" } grep { $aln[$_] =~/\bSS_cons\b/i } 1 .. @aln;
$ss_cons_len = @ss_cons;

#Check if the alignment is wrapped. It will be typically if the output is from mlocarna. This is important to check to calculate the motif 
my $motif_loc="";
my $ssConsensus=$ss_consLen="";
if($ss_cons_len == 1)
{
	my ($lineIdx, $string) = split(/\t/, $ss_cons[0]);
	$string =~ s/\s+/\t/g;
	($str1, $str2, $ssConsensus) = split(/\t/, $string);
	$ss_consLen = length($ssConsensus);
	
	#Get the motif start position in the alignment. By default we check for first open bracket: ( or { or [ or < - as an indication of a stem.
	$stem_open_brackets = qr/[\(\{\[<]/;
	$stem_close_brackets = qr/[\)\}\]>]/;
	
	my $stem_start_pos = index($ssConsensus, $1) if $ssConsensus =~ /($stem_open_brackets)/;
	my $stem_end_pos = rindex($ssConsensus, $1) if $ssConsensus =~ /($stem_close_brackets)/;

	$motif_loc = "$stem_start_pos\t$stem_end_pos";
}
else
{
	#Concatenate wrapped alignments in a line
	$ssConsensus="";
	@sorted = sort { $a <=> $b } @ss_cons; # Extra precautionary step to making sure lines are in the order, else it will goof up the alignment
	foreach $line (@sorted)
	{	
		my ($lineIdx, $string) = split(/\t/, $line);
		$string =~ s/\s+/\t/g;
		($str1, $str2, $p) = split(/\t/, $string);
		$ssConsensus = $ssConsensus.$p;
	}
	$ss_consLen = length($ssConsensus);
	
	#Get the motif start position in the alignment. By default we check for first open bracket: ( or { or [ or < - as an indication of a stem.
	$stem_open_brackets = qr/[\(\{\[<]/;
	$stem_close_brackets = qr/[\)\}\]>]/;
	
	my $stem_start_pos = index($ssConsensus, $1) if $ssConsensus =~ /($stem_open_brackets)/;
	my $stem_end_pos = rindex($ssConsensus, $1) if $ssConsensus =~ /($stem_close_brackets)/;

	$motif_loc = "$stem_start_pos\t$stem_end_pos";	
}

#Get the IDs of the sequences in the alignment
foreach $line (@aln)
{
        if($line !~ /^#|^\//)
        {
		$count +=1;
		next if($line !~ /[ATUGCatugc]/);
		$line =~ s/\s+/\t/g;
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
        foreach $key (@ids)
        {
                my $concat_seq = "";
                my @broken_seqs = grep { /\b$key\b/ } @aln_entries;
                my @sorted_brokenSeqs = sort { $a <=> $b } @broken_seqs;
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
        foreach $n (@aln_entires)
        {
                my ($srNo, $id, $p) = split(/\t/, $n);
                my $entry = $id."\t".$p;
                push(@alignments, $entry);
        }
}

#We take 10nt flanking from the motif location.
my ($motif_start, $motif_end) = split(/\t/, $motif_loc);

#If the motif starts at the start of the alignment. We wont have a 5' 10nt flank. Handling that with if..else
if($motif_start-10>1)
{
	$start = $motif_start-10;
}
else
{
	$start = 1;
}

#If the motif ends at the end of the alignment. We wont have a 3' 10nt flank. Handling that with if..else
if($motif_end+10<$ss_consLen)
{
	$end = $motif_end+10;
}
else
{
	$end = $ss_consLen;
}
my $motif_len = $end-$start+1;

#Construct output stockholm and fasta files for R-scape and RNAz respectively.
open(STO, ">$ARGV[1]"); #stockholm output file for R-scape
open(FA, ">$ARGV[2]"); #fasta output for RNAz
open(SQ, ">$ARGV[3]"); #fasta-like output for SQUARNA: For our purpose we are providing the locarna structure as reference with an alignment. Hence the format of this file is structure followed by aligned fasta sequences. Read SQUARNA manual for more details.

#Write the header from original results.stk file and for FA file
foreach $t (@top_lines)
{
	$t =~ s/^\s//;
	print STO $t;
}
print STO "#=GF Motif region with 10nt flank obtained with trimAlignment.pl\n\n";

#Get SS_cons from the motif region
my $ss_motifCons = substr($ssConsensus, $start, $motif_len);
print SQ "$ss_motifCons\n";

#Trim alignment
foreach my $entry (@alignments)
{
	#print "$entry\n\n";
	my ($id, $seq) = split(/\t/, $entry);
	my @values= split(/_/, $id);
	$win_chrInfo = join('_', @values[0 .. $#values - 3]);
	$win_end=$values[-1];
	$win_start=$values[-2];
	$win_counter=$values[-3];
		
	#Actual motif coordinates
	$new_start = $win_start+$start;
	$new_end = $win_start+$end;

	#Extract the motif from the alignment
	$m = substr($seq, $start, $motif_len);

	$identifier = "${win_chrInfo}_${win_counter}_${win_start}_${win_end}-${new_start}_${new_end}" if($win_end !~ /r\b/);
	$identifier = "${win_chrInfo}_${win_counter}_${win_start}_${win_end}-${new_start}_${new_end}r" if($win_end =~ /r\b/);
	print STO "$identifier\t\t$m\n";
	print FA ">$identifier\n$m\n";
	print SQ ">$identifier\n$m\n";
	undef($win_chrInfo), undef($win_start), undef($win_counter), undef($win_end);
}

#Printing SS_cons from the motif region
print STO "#=GC SS_cons\t\t$ss_motifCons\n";
print STO "//\n";

close(STO), close(FA), close(SQ);
