#!/usr/bin/perl 

###
# Pod Documentation
###

=head1 NAME

stockholm_to_html.pl

=head1 SYNOPSIS

Usage: stockholm_to_html.pl input.sto output.html

Create a colored HTML alignment from a Stockholm format sequence alignment.

=head1 DESCRIPTION

=over

=item B<--noLink>

Do not create hyperlinks from sequence accessions.

=back

=head1 AUTHOR

Zasha Weinberg
Jeffrey Barrick

=head1 COPYRIGHT

Copyright (C) 2008-2009 Zasha Weinberg, Jeffrey Barrick.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

###
# End Pod Documentation
###


use Bio::AlignIO;
use Getopt::Long;

$doPdf=0;
$rfamCsvFileName="";
$organismAssocFileName="";
$familyName="";
$maxHelixBulge=2; # >2 bases means we're not in a helix any more
%normalPairs=("AU",1,"UA",1,"CG",1,"GC",1,"GU",1,"UG",1);
$textwidthAdd=600;
$pageheight="6.5";
$ncbiSite = 1;

# first color is the one for bad base pairs
@colors=("#d0d0d0","#ff9999","#9999ff","#99ff99","#FF9900","#99ffff","#98C0C0","#ffff99","#ff33ff","#33ffff","#ffff33","#2ED70C","#F4BEF8","#ff9900","#B94F32","#FF0000","#ffcc99","#CCCCCC","#CC3366","#CCff66","#Ffcc66","#7e652f");

GetOptions( 
	    "h" => \$help,
	    "orgs=s" => \$organismAssocFileName,
	    "rfamCsv=s" => \$rfamCsvFileName,
	    "name=s" => \$familyName,
	    "pdf" => \$doPdf,
	    "nodup" => \$nodup,
	    "ncbiSite" => \$ncbiSite,
	    "textwidthAdd=s" => \$textwidthAdd,
	    "nonote" => \$nonote,
	    "pageheight=s" => \$pageheight,
	    "sortSpecies" => \$sortBySpecies,
	    "noLink" => \$noLink
	    );

if ($help) {
    print "perl stockholm_to_html.pl [options] <.sto file> <.html file> : create HTML file based on Stockholm-format file, with fancy formatting.\n";
    print "-orgs <file name> : file name is a tab-delim file that shows the relationship between embl Ids & organism.  Format is described in RaveNnA.pdf.  For EMBL sequecnes in their oroginal EMBL format, file can be made with ExtractEmblIdAndOrg2.pl\n";
    print "-rfamCsv <file name> : supply Rfam .csv file, which gives info on scores and on which seqs are new.\n";
    print "-name <text> : give the family a name; used in descriptive text.\n";
    print "-pdf : make PDF file instead of HTML file.  (Okay, in reality, it just makes a LaTeX file that you can compile with pdflatex.)\n";
    print "-nodup : don't print out duplicate sequences.  Sequences are duplicates iff they have the same organism name AND the same sequence.  Also, if sequences are the same, but one is NEW, both will be printed.  WARNING: currently only implemented for PDF.\n";
    print "-ncbiSite : make links go to NCBI site\n";
    print "-textwidthAdd <number of points> : add length in points to textwidth (only applicable to -pdf).  Default is 600 points.\n";
    print "-nonote : disable the note (description) on PDF, if you just want to print it out.\n";
    print "-pageheight <inches> : inches to set pdfpageheight, for pdflatex\n";
    print "-sortSpecies : sort the MSA by species name.  Currently this won't work if CM scores are supplied, since then the sort will be by score.  I'll work on it later maybe...\n";
    die;
}

$stoFileName=shift;
$htmlFileName=shift;

if ($stoFileName eq "" || $htmlFileName eq "") {
    die "didn't specify all params\n";
}

$explanation="Selected hits for $familyName are shown.  (Highly suspicious hits are removed, to avoid creating a confusing alignment.)  Hits that are not in Rfam are marked with the text \"NEW\" in bold at their left, and the nucleotides in their sequence are underlined.  The \"score\" column shows the CM-assigned score in bits.  The \"hit id\" column shows the location of the Hit in Rfam format: \"ID/start-end\", where ID is the EMBL nucleotide database ID of the sequence (which is hyperlinked to the EMBL database), and the nucleotides start..end are included in the hit sequence.  (If end is less than start, the hit is on the reverse strand.)  Nucleotides predicted to be in helices by the CM's Viterbi alignment are colored.  Base pairs in expected helices that are not Watson-Crick or G-U pairs are colored light gray; these may represent physical non-canonical base pairs, internal loops or a faulty alignment (all alignments were automatically generated, and only very clear imperfections were corrected).  Hits are listed in decreasing order of CM score.";
if ($nodup) {
    $explanation=$explanation."  NOTE: duplicate sequences in the same organism have been removed, so the actual number of hits is greater than the number shown.  More than one hit may have the sequence because an organism has more than one copy of the ncRNA, or simply because redundant sequences were entered into EMBL.  The duplicates were removed to make the alignment more readable regardless of whether the duplication has biological significance.\n";
}

use Data::Dumper;
#Open Input Alignment, now using BioPerl
my (@lines, @gc_lines, $ssSeq);
my $in = Bio::AlignIO->new(-format => 'stockholm', -file   => $stoFileName);
my $aln = $in->next_aln();
foreach my $seq ( $aln->each_seq() )
{
	push @lines, $seq->display_name . "/" . $seq->start . "-" . $seq->end . "_" . $seq->strand ."\t" . $seq->seq if($seq->strand eq "-1");
	push @lines, $seq->display_name . "/" . $seq->start . "-" . $seq->end ."\t" . $seq->seq if($seq->strand ne "-1");
}

my $ameta = $aln->consensus_meta;
foreach my $mname ($ameta->meta_names)
{
	if ($mname eq 'SS_cons')
	{
		$ssSeq = $ameta->named_meta($mname);
	}	
	push @gc_lines, "#=GC $mname " . $ameta->named_meta($mname);
}


## Older method of reading in -- does not handle interleaved format.
#open(STO,"$stoFileName");
# parse Stockholm file: read lines into @lines (tab delimit name\tseq), read $ssSeq
#while (<STO>) {
#    s/[\r\n]//g;
#
#    if ($_ eq "" || $_ eq "//") {
#	# ignore blank lines
#    }
#    else {
#
#	if (/^#/) {
#		if (/^#=GC/) {
#			push @gc_lines, $_;
#		    if (/^#=GC SS_cons/) {
#			/([^ ]*)$/;
#			$ssSeq=$1;
#			#print "$ssSeq\n";
#		    }
#	    }
#	}
#	else {
#	    /([^ ]*) [ ]*(.*)/;
#	    $s=uc $2;
#	    $line="$1\t$s";
#	    push @lines,$line;
#	}
#    }
#}

# for $ssSeq to use only angle brackets and periods, for my convenience
$_=$ssSeq;
s/\(/</g;
s/\)/>/g;
s/\[/</g;
s/\]/>/g;
s/\{/</g;
s/\}/>/g;
s/[-:,_]/./g;
#s/[a..z]/./g;
#s/[A..Z]/./g;

$ssSeq=$_;


my $done = 0;
my @pseudolines;
my $tempSeq = $ssSeq;
MAIN: while (!$done)
{
	for ($i=0; $i<length($tempSeq); $i++) {
	    $ch=substr $tempSeq,$i,1;

		if (!($ch =~ m/[<>.]/))
		{
			#print $ch ."\n";
			$upper = $ch;
			$lower = "\L$ch";
		
			#print "$lower$upper\n";
			my $char_throw = "[$lower$upper]";
			my $char_keep = "[^$lower$upper]";
			my $pseudo = $tempSeq;
			$tempSeq =~s/$char_throw/\./g;
			$pseudo =~ s/$char_keep/\./g;
			push @pseudolines, $pseudo;
		    #print "$pseudo\n$tempSeq\n";
			next MAIN;
		}
	}
	$done = 1;
}

# check for unexpected characters
#if (/[^.<>]/) {
#    die "Secondary structure line sequence contains unexpected characters.  After substitution of known characters, line is $ssSeq\n";
#}

# preprocess the $ssSeq -- for each position, map to the base-paired position
%ssPairMap=();
for ($i=0; $i<length($ssSeq); $i++) {
    $ch=substr $ssSeq,$i,1;

    $ssPairMap{$i}=-1;

    if ($ch eq "<") {
# look for matching partner
	$count=1;

	for ($j=$i+1; $j<length($ssSeq); $j++) {
	    $ch2=substr $ssSeq,$j,1;
	    if ($ch2 eq "<") {
		$count++;
	    }
	    if ($ch2 eq ">") {
		$count--;
	    }
	    if ($count==0) {
		$ssPairMap{$i}=$j;
		last;
	    }
	}

	if ($ssPairMap{$i}==-1) {
	    die "Secondary structure string doesn't match up.  $ssSeq\n";
	}
    }
}

#process pseudoknots
foreach my $pseudoline (@pseudolines)
{
	for ($i=0; $i<length($pseudoline); $i++) {
	    $ch=substr $pseudoline,$i,1;

	    if ($ch =~ m/[A-Z]/) {
	    
	 #   print "Char $ch\n";
	# look for matching partner
		$count=1;

		for ($j=$i+1; $j<length($pseudoline); $j++) {
		    $ch2=substr $pseudoline,$j,1;
		    if ($ch2 =~ m/[A-Z]/) {
			$count++;
		    }
		    if ($ch2 =~ m/[a-z]/) {
			$count--;
		    }
		   # print "$ch2 $count\n";
		    
		    if ($count==0) {
			$ssPairMap{$i}=$j;
			last;
		    }
		}

		if ($ssPairMap{$i}==-1) {
		    die "Secondary structure string doesn't match up.  $pseudoline\n";
		}
	    }
	}
}





# union with reverse map (mapping right to left)
%tempMap=%ssPairMap;
for $i (keys %tempMap) {
    if ($tempMap{$i}!=-1) {
	$ssPairMap{$tempMap{$i}}=$i;
    }
}

%hitIdToSeq=();
for $line (@lines) {
    ($hitId,$seq)=split /\t/,$line;
    $hitIdToSeq{$hitId}=$seq;
    print "$hitId\t$seq\n";
    push @hitIdList,$hitId;
}

#ignore mostly absent columns when giving helices color
my %not_mostly_gaps;
for (my $i=0; $i<length($ssSeq); $i++) {
	my $count = 0;
	for $hitId (@hitIdList) {
		$seq=$hitIdToSeq{$hitId};	

		my $c = substr $seq, $i, 1;
		$count++ if ($c =~ m/[ACTGUactgu]/);
	}
	
	$not_mostly_gaps{$i} = $count / (scalar @hitIdList) > 0.05;
}

#use Data::Dumper;
#print Dumper(\%not_mostly_gaps);


# for $i (keys %ssPairMap) { if ($ssPairMap{$i}!=-1) { print "$i --> $ssPairMap{$i}\n"; } }

# use heuristics to find helices & color them differently
%posToHelixNum=();
$currHelixNum=0;
for ($i=0; $i<length($ssSeq); $i++) {
    if ($ssPairMap{$i}!=-1 && $ssPairMap{$i}>$i) { # only do on left parts of helix
	$leftConnects=0;
	my $bulge_count = 0;
	for ($j=$i-1; $j>=0; $j--) {
		last if ($bulge_count>$maxHelixBulge);
		$bulge_count += $not_mostly_gaps{$j};
	    $ch=substr $ssSeq,$j,1;
	    if ($ch =~ m/[<A-Z]/) {
		$leftConnects=1;
	    }
	}
	$rightConnects=0;
	$bulge_count = 0;
	for ($j=$ssPairMap{$i}+1; $j<length($ssSeq); $j++) {
		last if ($bulge_count>$maxHelixBulge);
		$bulge_count += $not_mostly_gaps{$j};	
	    $ch=substr $ssSeq,$j,1;
	    if ($ch =~ m/[>a-z]/) {
		$rightConnects=1;
	    }
	}
	if ($leftConnects && $rightConnects) {
	    # same helix
	}
	else {
	    # new helix
	    $currHelixNum++;
	}

	$posToHelixNum{$i}=$currHelixNum;
	$posToHelixNum{$ssPairMap{$i}}=$currHelixNum;
    }
}

# parse the organisms, if requested
if (length($organismAssocFileName)>0) {
    if ($organismAssocFileName =~ /\.gz$/) {
	$organismAssocFileName="gzip -d -c $organismAssocFileName|";
    }
    open(ORGS,"$organismAssocFileName");

    #print "$organismAssocFileName\n";

    if (0) {
	# code for reading old format -- note: hasn't been tested with new compact internal representation
	$ocNum=1;
	while (<ORGS>) {
	    s/[\r\n]//g;
	    
	    ($emblId,$os,$oc)=split /\t/,$_;
	    
	    # just take first 2 words of species
	    $species=$os;
	    $_=$os;
	    if (/^([^ ]* [^ ]*)/) {
		$species=$1;
	    }
	    # represent this in the more compact internal format, but in an inefficient way
	    $ocNumToOc{$ocNum}=$species;
	    $emblIdToOcNum{$emblId}=$ocNum;
	    $ocNum++;
	}
    }
    else {
	# new format
	while (<ORGS>) {
	    s/[\r\n]//g;
	    if ($_ eq "") {
		$doingOc=1;
	    }
	    else {
		if ($doingOc) {
		    ($ocNum,$oc)=split /\t/,$_;
		    $ocNumToOc{$ocNum}=$oc;
		    #print "$ocNum,$oc\n";
		}
		else {
		    ($emblId,$ocNum)=split /\t/,$_;
		    $emblIdToOcNum{$emblId}=$ocNum;
		    #print "$emblId,$ocNum\n";
		}
	    }
	}
    }
    close(ORGS);
}



# parse RfamCsv
%hitIdToScore=();
%hitIdIsNew=();
if (length($rfamCsvFileName)>0) {
    open(CSV,"$rfamCsvFileName");
    $isFirst=1;
    while (<CSV>) {
	if ($isFirst) {
	    $isFirst=0;
	}
	else {
	    s/[\r\n]//g;
	    @fields=split /,/,$_;
	    $hitId=$fields[17];
	    $score=$fields[11];
	    $newness=$fields[15];
	    if ($newness==2) {
		$hitIdIsNew{$hitId}=1;
	    }
	    $hitIdToScore{$hitId}=$score;
	    push @scoreList,$score;
	    $scoreToHitIds{$score}=$scoreToHitIds{$score}."\t".$hitId;
	}
    }
    close(CSV);

    @scoreList=sort { $b <=> $a } @scoreList;
}
else {
    if ($sortBySpecies) {
	%scoreToHitIds=();
	@scoreList=("");

	@organismHitList=();
	for $hitId (@hitIdList) {
	    #print "Q: $hitId\n";
	    $hitId =~ /([-._A-Za-z0-9]*)\/[0-9]*[-][0-9]*/;
	    $emblId=$1;
	    $species=EmblIdToSpecies($emblId);
	    $t="$species\t$hitId";
	    #print "$hitId,$emblId,$species,$t\n";
	    push @organismHitList,$t;
	}
	@organismHitList=sort { (lc $a) cmp (lc $b) } @organismHitList;
	$t="";
	for $orgHit (@organismHitList) {
	    $orgHit =~ /([^\t]*)\t(.*)/;
	    $org=$1;
	    $hitId=$2;
	    $t=$t."\t".$hitId;
	}
	#print "WWW: $t\n";
	$scoreToHitIds{""}=$t;
    }
    else {
	%scoreToHitIds=();
	$t="";
	for $hitId (@hitIdList) {
	    $t=$t."\t".$hitId;
	}
	$scoreToHitIds{""}=$t;
	@scoreList=("");
    }
}

# dump the HTML or PDF

if ($doPdf) {
    
# create PDF via LaTeX

    if ($htmlFileName =~ /\.tex$/) {
	$texFileName=$htmlFileName;
    }
    else {
	$texFileName="$htmlFileName.tex";
    }
    open(LATEX,">$texFileName");
    print LATEX "\\documentclass{article}\n\\usepackage{color}\n\\setlength{\\fboxsep}{0pt}\n";
    print LATEX "\\usepackage{longtable}\n";
    print LATEX "\\usepackage{amsmath}\n";
    
    $i=0;
    for $color (@colors) {
	@component=();
	push @component,lc(substr $color,1,2);
	push @component,lc(substr $color,3,2);
	push @component,lc(substr $color,5,2);
	# print "$color,@component\n";

	print LATEX "\\definecolor{b$i}{rgb}{";
	for ($j=0; $j<3; $j++) {
	    if ($j>0) {
		print LATEX ",";
	    }
	    $level=hex($component[$j])/255.0;
	    print LATEX "$level";
	}
	print LATEX "}\n";
	$i++;
    }
    
    print LATEX "%pdflatex params & extraction from hyperref package\n";
    print LATEX "\\pdfhorigin0pt\n\\pdfvorigin0pt\n\\pdfpagewidth0pt\n";
    print LATEX "\\pdfpageheight".$pageheight."in\n";
    print LATEX "\\newcommand{\\myhref}[2]{\\pdfstartlink attr{/C [0 0 0.9] /Border [0 0 1]} user{/Subtype/Link/A<</Type/Action/S/URI/URI(#1)>>} #2 \\pdfendlink }\n";

    print LATEX "\\addtolength\\textwidth{$textwidthAdd"."pt}\n\\addtolength\\textwidth{72pt}\n";
    $textheight=$pageheight-2;
    print LATEX "\\setlength\\textheight{$textheight"."in}\n";
    print LATEX "\\pagestyle{empty}\n";

    print LATEX "\\begin{document}\n";

    if (!$nonote) {
	$_=$explanation;
	s/\"([^\"]*)\"/``$1''/g;
	s/underlined/boxed/g;
	print LATEX "\\parbox{7in}{$_}\\\\\n";
    }

    print LATEX "\\vspace*{3ex}\n";

    print LATEX "\\begin{longtable}{llll}\n";

    # secondary structure string goes at the bottom of each page
    $_=$ssSeq;
    s/</&lt;/g;
    s/>/&gt;/g;
    print LATEX "&&&{\\tt $ssSeq}\\\\\n";
    print LATEX "\\endfoot\n";

    # headers
    if (length($rfamCsvFileName)>0) {
	print LATEX "{\\bf Score (bits)}&{\\bf Sequence/Start-End}&";
    }
    else {
	print LATEX "&&";
    }
    if (length($organismAssocFileName)>0) {
	print LATEX "{\\bf Species}&";
    }
    else {
	print LATEX "&";
    }
    print LATEX "{\\bf Sequence}\\\\\n";

    %seqAndOrgSeen=(); # do avoid duplicates, if desired

    $prevScore="qqqqqq"; # impossible score string
    for $score (@scoreList) {

	#print "$score,$prevScore\n";

	if ($score eq $prevScore) {
	    next;
	}
	$prevScore=$score;

	$hitIds=$scoreToHitIds{$score};
	@hitIdList=split /\t/,$hitIds;
	for $hitId (@hitIdList) {
	    if (length($hitIdToSeq{$hitId})>0) {
		$seq=$hitIdToSeq{$hitId};
		    
		$_=$hitId;
		$emblId="";
		$start="";
		$end="";
		if (/([.A-Za-z0-9_]*)\/([0-9]*)[-]([0-9]*)/) {
		    $emblId=$1;
		    $start=$2;
		    $end=$3;
		}
		
		$showThisLine=1;

		# check for duplicates, if desired
		if ($nodup) {
		    $species=EmblIdToSpecies($emblId);
		    $dupKey="$hitIdIsNew{$hitId}\t$species\t$seq";
		    if ($seqAndOrgSeen{$dupKey}) {
			$showThisLine=0;
		    }
		    else {
			$seqAndOrgSeen{$dupKey}=1;
		    }
		    #print "$score, $dupKey: $showThisLine\n";
		    #print %seqAndOrgSeen;
		    #print "\n\n";
		}
		
		if ($showThisLine) {

		    if ($hitIdIsNew{$hitId}) {
			print LATEX "{\\bf NEW}\\hfill ";
		    }
		    
		    print LATEX "$score&";
		    
		    if ($emblId eq "") {
			print LATEX "$hitId&";
		    }
		    else {
			$_=$emblId;
			/(^[^.]*)/;
			$_=$1;
			s/[_]/\\_/g;
			$strippedEmblId=$_;
			$_=$emblId;
			s/_/\\_/g;
			$escapedEmblId=$_;
			if ($ncbiSite) {
			    print LATEX "\\myhref{http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=nucleotide&orig_db=nucleotide&term=$strippedEmblId}{$escapedEmblId}/$start-$end&";
			}
			else {
			    print LATEX "\\myhref{http://www.ebi.ac.uk/cgi-bin/emblfetch?style=html&id=$strippedEmblId}{$escapedEmblId}/$start-$end&";
			}
		    }
		    
		    $species=EmblIdToSpecies($emblId);
		    print LATEX "$species&";
		    
		    print LATEX "{\\tt ";
		    if ($hitIdIsNew{$hitId}) {
			print LATEX "\$\\underline{\\text{";
		    }
		    
		    for ($i=0; $i<length($seq); $i++) {
			$ch=substr $seq,$i,1;
			
			if ($ssPairMap{$i}==-1) {
			    print LATEX "$ch";
			}
			else {
			    if ($ch eq "-") {
				# don't highlight deletants so much
				print LATEX "$ch";
			    }
			    else {
				$ch2=substr $seq,$ssPairMap{$i},1;
				$pair="$ch$ch2";
				if ($normalPairs{$pair}) {
				    $helixNum=$posToHelixNum{$i};
				    print LATEX "\\colorbox{b$helixNum}{$ch}";
				}
				else {
				    print LATEX "\\colorbox{b0}{$ch}";
				}
			    }
			}
		    }
		    if ($hitIdIsNew{$hitId}) {
			print LATEX "}}\$";
		    }
		    
		    print LATEX "}\\\\\n";
		}
	    }
	}
    }

    print LATEX "\\end{longtable}\n";
    print LATEX "\\end{document}\n";
}
else {
# HTML
    open(HTML,">$htmlFileName");
    
    print HTML "<html><head>";
    if ($familyName eq "") {
    }
    else {
	print HTML "<title>$familyName: multiple alignment highlighting new hits with rigorous filtering</title>";
    }
    if (1) {
# inline style
	print HTML "<style>b, .b { font-weight: normal;}\n";
	print HTML "#reallyBold { font-weight:bold; }\n";
    $i=0;
    for $color (@colors) {
	print HTML "#b$i {background-color: $color;}\n";
	$i++;
    }
    print HTML "#bad {background-color: $color[0];}\n";
	print HTML "</style>\n";
    }
    else {
	# .css file
	print HTML "<link rel=stylesheet href=\"./FancifyStockholm.css\" type=\"text/css\">";
    }
    print HTML "</head><body>\n";
    if ($familyName eq "") {
    }
    else {
	# this comment is somewhat specific to the NAR paper
	$_=$explanation;
	s/\"/&quot;/g;
	print HTML "<p>$_</p>\n";
    }
    print HTML "<table cellspacing=0 cellpadding=0 border=0>\n";
    print HTML "<tr><td></td>";
    if (length($rfamCsvFileName)>0) {
	print HTML "<td><b ID=\"reallyBold\">Score (bits)</b></td>";
    }
    else {
	print HTML "<td></td>";
    }
    print HTML "<td><b ID=\"reallyBold\">Accession/Start-End</b></td><td>";
    if (length($organismAssocFileName)>0) {
	print HTML "<b ID=\"reallyBold\">Species</b>";
    }
    print HTML "</td><td><b ID=\"reallyBold\">Sequence</b></td>";
    print HTML "</tr>\n";

    $prevScore="qqqqqq"; # impossible
    for $score (@scoreList) {

	if ($score eq $prevScore) {
	    next;
	}
	$prevScore=$score;

	$hitIds=$scoreToHitIds{$score};
	@hitIdList=split /\t/,$hitIds;
	for $hitId (@hitIdList) {
	    if (length($hitIdToSeq{$hitId})>0) {
		$seq=$hitIdToSeq{$hitId};
		
		$_=$hitId;
		$emblId="";
		$start="";
		$end="";
		if (/([.A-Za-z0-9_]*)\/([0-9]*)[-]([0-9]*)/) {
		    $emblId=$1;
		    $start=$2;
		    $end=$3;
		}
		
		$species=EmblIdToSpecies($emblId);
		$onMouseOver="JavaScript:window.status='$hitId    $species'";
		
		print HTML "<tr>";
		if ($hitIdIsNew{$hitId}) {
		    print HTML "<td><b ID=\"reallyBold\">NEW</b></td>";
		}
		else {
		    print HTML "<td></td>";
		}
		print HTML "<td><nobr>$score</nobr></td>";
		if ($emblId eq "") {
		    print HTML "<td><nobr>$hitId</nobr></td>";
		}
		else {
		    $_=$emblId;
		    /(^[^.]*)/;
		    $strippedEmblId=$1;
		    if ($noLink)
		    {
		   	 	print HTML "<td><nobr>$emblId/$start-$end</nobr>";
		   	}
		   	else
		    {
		    	if ($ncbiSite) {
				print HTML "<td><nobr><a href=\"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=nucleotide&orig_db=nucleotide&term=$strippedEmblId\">$emblId</a>/$start-$end</nobr>&nbsp;&nbsp;</td>";
			    }
			    else {
				print HTML "<td><nobr><a href=\"http://www.ebi.ac.uk/cgi-bin/emblfetch?style=html&id=$strippedEmblId\">$emblId</a>/$start-$end</nobr>&nbsp;&nbsp;</td>";
			    }
			}
		}
		$species=EmblIdToSpecies($emblId);
		print HTML "<td><nobr>$species</nobr></td>";
		print HTML "<td onMouseOver=\"$onMouseOver\" onMouseOut=\"JavaScript:window.status=''\"><tt><nobr>";
		if ($hitIdIsNew{$hitId}) {
		    print HTML "<u>";
		}

		for ($i=0; $i<length($seq); $i++) {
		    $ch=substr $seq,$i,1;
		    
		    if ($ssPairMap{$i}==-1) {
			print HTML "$ch";
		    }
		    else {
			if ($ch eq "-") {
			    # don't highlight deletants so much
			    print HTML "$ch";
			}
			else {
			    $ch2=substr $seq,$ssPairMap{$i},1;
			    $pair="$ch$ch2";
			    if ($normalPairs{$pair}) {
				$helixNum=$posToHelixNum{$i};
				print HTML "<b id=\"b$helixNum\">$ch</b>";
			    }
			    else {
				print HTML "<b id=\"bad\">$ch</b>";
			    }
			}
		    }
		}
		if ($hitIdIsNew{$hitId}) {
		    print HTML "</u>";
		}

		print HTML "</nobr></tt>";
		print HTML "</td></tr>\n";
	    }
	}
    }

    $_=$ssSeq;
    s/</&lt;/g;
    s/>/&gt;/g;
    
    foreach my $line (@gc_lines)
    {
    	$line =~ m/#=GC\s+(\S+)\s+(\S+)/;
    	
    	my $title = $1;
    	my $content = $2;
    	$content =~ s/</&lt;/g;
    	$content =~ s/>/&gt;/g;

    	print HTML "<tr><td></td><td></td><td>$title</td><td></td><td><tt>$content</tt></td></tr>\n";
    }
    print HTML "</table></body></html>\n";
}

sub EmblIdToSpecies {
    my ($emblId)=@_;
    my $ocNum=$emblIdToOcNum{$emblId};
    my $oc=$ocNumToOc{$ocNum};
    my @ocList=split /;/,$oc;
    my $species=pop @ocList;
    #print "EmblIdToSpecies: $emblId,$ocNum,$oc,@ocList,$species\n";
    return $species;
}
