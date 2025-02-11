#! usr/bin/perl
use POSIX;
use Getopt::Long;
use Pod::Usage;

GetOptions(
        'g|genomeid=s'	=> \$id,	#Genome ID (Typically NC ID from NCBI)
        's|scores=s'	=> \$sfile,	#trailing info in the filename (Example: if filename is NC_007407.1_102_17850_18100r_scores.tab, the -s option should be given as -s _scores.tab
        'e|eval=s'	=> \$efile,	#trailing info in the filename (Example: if filename is NC_007407.1_102_17850_18100r_hits_evaluated.tsv, the -s option should be given as -s _hits_evaluated.tsv
	'h|help'	=> \$help,
	) or die "Usage: Incorrect Usage.\n Use --help to see the correct use of the script\n";

#Display help if needed
pod2usage(1) if $help;

# Check if the mandatory extension option is provided
if(!defined $ARGV[0])
{
        die "Error: Input files are not given\nProvide the script with -g option and 1 input file\nExample usage: perl cutsv2bed.pl -g NC_000932 -s _scores.tab -e _hits_evaluated.tsv input_clusterList.csv \n";
}
if (!defined $id or !defined $sfile or !defined $efile )
{
        die "Error: The genomeID (-g option) and trailing file extensions (-s and -e options) are mandatory.\nExample usage: perl cutsv2bed.pl -g NC_000932 -s _scores.tab -e _hits_evaluated.tsv input_clusterList.csv \n";
}

#Reading cluster list file
open(CLUST, $ARGV[0]);
my @clust = <CLUST>;
close(CLUST);

print "track name=\"RNAmotifs\" description=\"RNA motifs identified through extensive computational analyses\" visibility=2 colorByStrand=\"200,90,222 238,187,250\" ggfTags=on\n"; 
foreach $clust (@clust)
{
	chomp($clust);

	#Reading scores and evaluation file 
	$score_file=${clust}."$sfile"; #Modify according to you filenames. Note the filenames containing more than cluster ID (from cluster file).
	open(SCORES, $score_file);
	my @scores = <SCORES>;
	close(SCORES);

	$eval_file=${clust}."$efile"; #Modify according to you filenames. Note the filenames containing more than cluster ID (from cluster file).
	open(EVAL, $eval_file);
	my @eval = <EVAL>;
	close(EVAL);

	#Getting homolog information for genome ID from the evaluation file
	foreach $eval_entry (@eval)
	{
		my $seqRep, $length, $GCcontent, $bp_threshold, $total_bp, $openBP, $closeBP, $hairpinPER, $bpActuallyFormed, $RedundantRNAs, $evaluation, @match;
		my $genomeID, $RNAcoords, $RNAstart, $RNAend, $RNAstrand, $desc, $bitscore, $evalue;
		if($eval_entry =~ /$id/ && $eval_entry =~ /Passed/)
		{
			($seqRep, $length, $GCcontent, $bp_threshold, $total_bp, $openBP, $closeBP, $hairpinPER, $bpActuallyFormed, $RedundantRNAs, $evaluation) = split(/\t/, $eval_entry);
			if($seqRep =~ /$id/)
			{
				my $match = $seqRep;
				push(@match, $match);
			}
			if($RedundantRNAs =~ /$id/)
			{
				my $match = $RedundantRNAs;
				foreach my $RNA_identHomolog (split /,/, $match)
				{
					if($RNA_identHomolog =~ /$id/)
					{
						push(@match, $RNA_identHomolog);
					}
				}
			}	

			foreach my $entry_match (@match)
			{
				($genomeID, $RNAinfo) = split(/\//, $entry_match);
				($RNAcoords, $RNAstrand) = split(/_/, $RNAinfo);
				($RNAstart, $RNAend) = split(/-/, $RNAcoords);
	
				#Get motif scores
				foreach $score (@scores)
				{	
					my $clusterID, $RNA, $strandInfo;
					my $start, $end;
					if($RNAstart > $RNAend)
					{
						$start = $RNAend;
						$end = $RNAstart;
					}
					else
					{
						$start = $RNAstart;
						$end = $RNAend;
					}	
					
					if($score =~ /${genomeID}\/${start}-${end}/)
					{
						($clusterID, $desc, $RNA, $strandInfo, $bitscore, $evalue) = split(/\t/, $score);
						$desc =~ s/chloroplast, complete genome//g;
						chomp($desc);
					}
				}
				chomp($evalue);

				#Generate BED file
				$strand = $RNAstrand;
				$RNAstrand = ($RNAstrand == 1) ? "+" : "-";
				$strand = ($strand == 1) ? "plus" : "minus";
				$name="RNAcoords=$genomeID/$RNAstart-$RNAend,RNAstrand=$strand,"."SequenceRep=$seqRep,"."model=$clust,"."bit_score=$bitscore,"."E-value=$evalue,"."Organism=$desc";
				chomp($name);
				print "$genomeID\t$RNAstart\t$RNAend\t$name\t0\t$RNAstrand\t$RNAstart\t$RNAend\n";
			}
			undef(@match);
		}
	}

	undef(@scores), undef(@eval);
}

undef(@clust);
