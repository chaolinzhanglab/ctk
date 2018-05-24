#!/usr/bin/perl -w

use strict;
use warnings;

use File::Basename;
use Getopt::Long;

use MyConfig;
use Bed;
use ScanStatistics;

use Data::Dumper;

my $prog = basename ($0);
my $cmdDir = dirname ($0);

my $verbose = 0;
my $big = 0;
my $minBlockSize = 2000000;
my $separateStrand = 0;

my $valleySeeking = 0;	#0 or 1
my $valleyDepth = 0.9;	#0.9 as default

#my $confFile = "";
my $dbkey = "";

my $geneBedFile = "";
my $useExpr = 0;
my $pvalueThreshold = 0.01;
my $multiTestCorrection = 0;

my $minPH = 2;
my $maxPH = -1;
my $skipOutofRangePeaks = 0; #by default, peaks with PH> maxPH will be labeled with p=nd, but kept. This option is provided to reproduce results from a previous version

#calculation of scan statistics can be slow for very large peak height, and these might be spurious (e.g., in very abundant RNA) and one might filter them out
#5000; -1= no filter by default

my $maxGap = -1; #no merge by default, default 25 recommended if merge

my $outBoundaryBedFile = "";
my $outHalfPHBedFile = "";

my $prefix = "Peak";

my $cache = getDefaultCache ($prog);
my $keepCache = 0;

my @ARGV0 = @ARGV;

GetOptions (
	'big'=>\$big,
	'p:f'=>\$pvalueThreshold,
	'ss'=>\$separateStrand,
	'valley-seeking'=>\$valleySeeking,
	'valley-depth:f'=>\$valleyDepth,
	#'test'=>\$doStatTest,
	'dbkey:s'=>\$dbkey,
	'gene:s'=>\$geneBedFile,
	'use-expr'=>\$useExpr,
	'multi-test'=>\$multiTestCorrection,
	'minPH:i'=>\$minPH,
	'maxPH:i'=>\$maxPH,
	'skip-out-of-range-peaks'=>\$skipOutofRangePeaks,
	'gap:i'=>\$maxGap,
	'out-boundary:s'=>\$outBoundaryBedFile,
	'out-half-PH:s'=>\$outHalfPHBedFile,
	'prefix:s'=>\$prefix,
	'c|cache:s'=>\$cache,
	'keep-cache'=>\$keepCache,
	'v'=>\$verbose);


if (@ARGV != 2)
{
    print "detecting peaks from CLIP data\n";
    print "Usage: $prog [options] <tag.bed> <peak.bed>\n";
    print " <tag.bed> : BED file of unique CLIP tags, input\n"; #, use - for stdin\n";
	print " <peak.bed>: BED file of called peaks, output\n";
	print "Options:\n";
	#print " -e          : expression values indicated in the gene bed file\n";
    print " -big                   : big input file\n";
	print " -ss                    : separate the two strands\n";
	print " --valley-seeking       : find candidate peaks by valley seeking\n";
	print " --valley-depth [float] : depth of valley if valley seeking ($valleyDepth)\n";
	print " --out-boundary [string]: output cluster boundaries\n";
	print " --out-half-PH  [string]: output half peak height boundaries\n";
    print " --dbkey        [string]: species to retrieve the default gene bed file (mm10|hg19)\n";
	print " --gene         [string]: custom gene bed file for scan statistics (will override --dbkey)\n";
	print " --use-expr             : use expression levels given in the score column in the custom gene bed file for normalization\n";
	print " -p             [float] : threshold of p-value to call peak ($pvalueThreshold)\n";    
	print " --multi-test           : do Bonferroni multiple test correction\n";
	print " -minPH         [int]   : min peak height ($minPH)\n";
    print " -maxPH         [int]   : max peak height to calculate p-value($maxPH, no limit if < 0)\n";
	print " --skip-out-of-range-peaks: skip peaks with PH > maxPH\n";
	print " -gap           [int]   : merge cluster peaks closer than the gap ($maxGap, no merge if < 0)\n";
	print " --prefix       [string]: prefix of peak id ($prefix)\n";
    print " -c             [dir]   : cache dir\n";
	print " --keep-cache           : keep cache when the job done\n";
    print " -v                     : verbose\n";
    exit (0);
}

print "CMD=$prog ", join(" ", @ARGV0), "\n" if $verbose;

my ($tagBedFile, $outBedFile) = @ARGV;


#check parameters

if ($geneBedFile ne '' || $useExpr)
{
	Carp::croak "$geneBedFile does not exist\n" unless -f $geneBedFile;
}


#the lines below is expired because we check this later
#if ($dbkey ne '')
#{
#	Carp::croak "dbkey must be mm10 or hg19\n" unless $dbkey eq 'mm10' || $dbkey eq 'hg19';
#}


if ($valleySeeking)
{
	Carp::croak "valley-depth must be between 0.5 and 1\n" unless $valleyDepth >= 0.5 && $valleyDepth <= 1;
}
else
{
	#if no valley seeking, must specify genic region and do statistical test
	#also must do merging
	Carp::croak "must specify --dbkey or --gene when --valley-seeking is not unabled\n" unless $geneBedFile ne '' || $dbkey ne '';
	warn "peak merging is strongly recommended (e.g. --gap 25) without valley seeking\n" if $maxGap < 0;
}

if ($dbkey ne '')
{
	if ($geneBedFile eq '')
	{
		#get gene bed file
		my $confFile = "$cmdDir/ctk.loc";
		my $locationInfo = getLocationInfo ($confFile, $dbkey, "genic");
		Carp::croak "cannot locate genic bed file for $dbkey\n" unless exists $locationInfo->{'genic'};

		$geneBedFile = $locationInfo->{'genic'};
		print "gene bed file=$geneBedFile\n" if $verbose;

		Carp::croak "$geneBedFile does not exist\n" unless -f $geneBedFile;
	}
	#otherwise do nothing
}

if ($outBoundaryBedFile ne '' || $outHalfPHBedFile ne '')
{
	Carp::croak "must enable --valley-seeking\n" unless $valleySeeking;
}

if ($minPH < 2)
{
	warn "minPH must be >=2\n";
	$minPH = 2;
}


my $bigFlag = $big ? "-big" : "";
my $verboseFlag = $verbose ? "-v" : "";
my $ssFlag = $separateStrand ? "-ss" : "";

Carp::croak "$cache already exists\n" if -d $cache || -f $cache;
system ("mkdir $cache");


##
print "generate exact tag coverage profile ...\n" if $verbose;

my $tagExactCountBedFile = "$cache/tag.count.exact.bed";
my $tmpPeakBedFile = "$cache/tmp.peak.bed";

my $cmd = "perl $cmdDir/tag2profile.pl $verboseFlag $bigFlag -exact $ssFlag -of bed $tagBedFile $tagExactCountBedFile";
my $ret = system ($cmd);
Carp::croak "CMD=$cmd failed:$?\n" unless $ret == 0;


print "extract candidate peaks ...\n" if $verbose;

if ($valleySeeking)
{
	my %tagCount;
	my $cachePeak = "$cache/cache_peak";
	system ("mkdir $cachePeak");
	if ($big)
	{
		my $ret = splitBedFileByChrom ($tagExactCountBedFile, $cachePeak, v=>$verbose, "sort"=>1);
		%tagCount = %$ret;
	} 
	else
	{
		my $tags = readBedFile ($tagExactCountBedFile, $verbose);
	   	foreach my $t (@$tags)
    	{
        	my $chrom = $t->{"chrom"};
        	push @{$tagCount{$chrom}}, $t;
    	}
	}

	print "get positions with non-zero tag count broken down into chromosomes ...\n" if $verbose;

	foreach my $chrom (sort keys %tagCount)
	{

	    my $n = $tagCount{$chrom};
	    $n = ref($n) eq 'HASH' ? $n = $n->{'n'} : @$n;
	    print "$chrom : $n tags\n" if $verbose;
	}

	my @strand = ('b');
	@strand = qw(+ -) if $separateStrand;

	my $fout;
	open ($fout, ">$tmpPeakBedFile") || Carp::croak "cannot open file $tmpPeakBedFile to write\n";
	foreach my $s (@strand)
	{
		print "#processing strand $s ...\n" if $verbose;
		my @chroms = sort keys %tagCount;

		foreach my $chrom (@chroms)
		{
			print "processing chrom $chrom ...\n" if $verbose;
			if ($big)
			{
				my $tmpFile = $tagCount{$chrom}->{'f'};
				my $fin;
				open ($fin, "<$tmpFile") || Carp::croak "cannot open file $tmpFile to read\n";

				my $iter = 0;
				while (my $block = readNextBedBlock ($fin, 1, 0, minBlockSize=> $minBlockSize))
				{
					my $n = @$block;
					print "block $iter: $n positions\n" if $verbose;
					$iter++;
					my $peaks = extractPeaks ($block, $s, $valleyDepth);
					foreach my $p (@$peaks)
					{
						$p->{'name'} .= "#" . join ("#", $p->{'leftBoundary'}, $p->{'rightBoundary'}+1, $p->{'leftHalfPH'}, $p->{'rightHalfPH'}+1);
						print $fout bedToLine ($p), "\n" if $p->{'score'} >= $minPH;
					}
				}
				close ($fin);
			}
			else
			{
				my $tags = $tagCount{$chrom};
				my $n = @$tags;
				print "$n tags loaded on chromosome $chrom\n" if $verbose;
				my $peaks = extractPeaks ($tags, $s, $valleyDepth);

   	        	foreach my $p (@$peaks)
            	{
					$p->{'name'} .= "#" . join ("#", $p->{'leftBoundary'}, $p->{'rightBoundary'}+1, $p->{'leftHalfPH'}, $p->{'rightHalfPH'}+1);
            		print $fout bedToLine ($p), "\n" if $p->{'score'} >= $minPH;
            	}
			}
		}
	}
	close ($fout);
}
else
{
	my $cmd = "awk '{if (\$5>=$minPH) {print \$0}}' $tagExactCountBedFile > $tmpPeakBedFile";
	system ($cmd);
}


$cmd = "wc -l $tmpPeakBedFile | awk '{print \$1}'";
my $n = `$cmd`;
$n=~/^(\d+)/;
$n = $1;


if ($n == 0)
{
	Carp::croak "no peaks found\n";
	system ("rm -rf $cache") unless $keepCache;
	exit (0);
}

print "$n candidate peaks with PH>=$minPH detected\n" if $verbose;

#assign candidate peaks to genes and assess statistical significance

if (-f $geneBedFile)
{
	print "estimating average tag size ...\n" if $verbose;
	my $tagSize =  `awk 'BEGIN{s=0;n=0;} {s=s+\$3-\$2; n=n+1} END {print s/n}' $tagBedFile`;
	
	chomp $tagSize;
	print "average tag size = $tagSize\n" if $verbose;

	print "get exons in genes ... \n" if $verbose;
	my $ts2geneFile = "$cache/ts2gene.txt";
	my $cleanGeneBedFile = "$cache/gene.clean.bed";

	my $cmd = "awk '{print \$4\"\\t\"\$4}' $geneBedFile | sort | uniq > $ts2geneFile";
	system ($cmd);

	$cmd = "perl $cmdDir/combineTranscripts.pl $verboseFlag  $geneBedFile $ts2geneFile $cleanGeneBedFile";
	system ($cmd);

	my $exonBedFile = "$cache/exon.bed";
	$cmd = "perl $cmdDir/gene2ExonIntron.pl $verboseFlag -nid -oe $exonBedFile $cleanGeneBedFile";
	system ($cmd);


	print "count tag number for each exon/gene ...\n" if $verbose;
	my $exonTagCountBedFile = "$cache/tag.count.exon.bed";

	$cmd = "perl $cmdDir/tag2profile.pl $bigFlag $verboseFlag -region $exonBedFile $ssFlag -of bed $tagBedFile $exonTagCountBedFile";
	print $cmd, "\n";
	system ($cmd);

	my $fin;
	my $fout;

	my %geneTagCountHash;

	open ($fin, "<$exonTagCountBedFile") || Carp::croak "cannot open file $exonTagCountBedFile\n";
	while (my $line = <$fin>)
	{
	    chomp $line;
	    my $r = lineToBed ($line);

	    my $geneId = $r->{"name"};
	    $geneTagCountHash{$geneId}->{'size'} += ($r->{'chromEnd'} - $r->{'chromStart'} + 1);
	    $geneTagCountHash{$geneId}->{'count'} += $r->{'score'};
	}
	close ($fin);

	my $geneTagCountTotal = 0;
	map {$geneTagCountTotal += $geneTagCountHash{$_}->{'count'}} keys %geneTagCountHash;

	Carp::croak "no tags??" if $geneTagCountTotal <= 0;

	print "total number of tags overlapping with specified reions: $geneTagCountTotal\n" if $verbose;

	my $geneExpressionTotal = 0;
	if ($useExpr)
	{
		print "reading gene expression levels from $geneBedFile ...\n" if $verbose;
		open ($fin, "<$geneBedFile") || Carp::croak "cannot open file $geneBedFile to read\n";
		while (my $line = <$fin>)
		{
			chomp $line;
			next if $line =~/^\s*$/;
			next if $line =~/^track/;
			next if $line =~/^#/;

			my @cols = split (/\s+/, $line);
			Carp::croak "no score column in $geneBedFile ...\n" if @cols < 5;
		
			my $geneId = $cols[3];
			my $expr = $cols[4];
		
			Carp::croak "expression value for gene $geneId is negative\n" if $expr < 0;
			$geneTagCountHash{$geneId}->{'expr'} = $expr if exists $geneTagCountHash{$geneId};
		}
		close ($fin);

		map {$geneExpressionTotal += $geneTagCountHash{$_}->{'expr'}} keys %geneTagCountHash;
		Carp::croak "all genes have no expression??\n" if $geneExpressionTotal <= 0;

		print "total gene expression level: $geneExpressionTotal\n" if $verbose;
	}

	print "calculating expected peak height for each gene ...\n" if $verbose;


#my $geneTagCountFile = "$cache/tag.count.gene.txt";
#open ($fout, ">$geneTagCountFile") || Carp::croak "cannot open file $geneTagCountFile to write\n";

#if ($useExpr)
#{
#	print $fout "#", join ("\t", "gene id", "size", "tag number", "tag density", "expression", "expexpted PH"), "\n";
#}
#else
#{
#	print $fout "#", join ("\t", "gene id", "size", "tag number", "tag density", "expexpted PH"), "\n";
#}

	foreach my $geneId (sort keys %geneTagCountHash)
	{
		my $g = $geneTagCountHash{$geneId};
	    $g->{'density'} = $g->{'count'} / $g->{'size'};
		$g->{'PH0'} = $useExpr ? ($geneTagCountTotal * $g->{'expr'} / $geneExpressionTotal / $g->{'size'} * $tagSize) : ($g->{'density'} * $tagSize);

#	if ($useExpr)
#	{
#		print $fout join ("\t", $geneId, $g->{'size'}, $g->{'count'}, $g->{'density'}, $g->{'expr'}, $g->{'PH0'}), "\n";
#	}
#	else
#	{
#		print $fout join ("\t", $geneId, $g->{'size'}, $g->{'count'}, $g->{'density'}, $g->{'PH0'}), "\n";
#	}

	}
#close ($fout);

#my $effectiveGeneNum = `awk '{if(\$3>0) {print \$0}}' $geneTagCountFile | wc -l | awk '{print \$1}'`;
#$effectiveGeneNum =~/\s*(\d+)\s*/;
#$effectiveGeneNum = $1;

	my $effectiveGeneNum = 0;

	if ($useExpr)
	{
		map {$effectiveGeneNum++ if $geneTagCountHash{$_}->{'expr'} > 0} keys %geneTagCountHash;
	}
	else
	{
		map {$effectiveGeneNum++ if $geneTagCountHash{$_}->{'count'} > 0} keys %geneTagCountHash;
	}

	print "effective gene number: $effectiveGeneNum\n" if $verbose;


	#this part can be improved
	#e.g. 1. if two genes overlap with each other resulting in ambiguity, exonic region should have higher priority
	#     2. peaks in the extended 3' UTR region

	print "match peaks with genes ...\n" if $verbose;
	my $peak2geneBedFile = "$cache/peak2gene.bed";
	$cmd = "perl $cmdDir/tagoverlap.pl $bigFlag $verboseFlag -d \"@@\" --keep-score -region $cleanGeneBedFile $ssFlag $tmpPeakBedFile $peak2geneBedFile";
	system ($cmd);


	print "calculate pvalues ...\n" if $verbose;

	$tmpPeakBedFile = "$cache/tmp.peak.sig.bed";
	open ($fin, "<$peak2geneBedFile") || Carp::croak "cannot open file $peak2geneBedFile to read\n";
	open ($fout, ">$tmpPeakBedFile") || Carp::croak "cannot open file $tmpPeakBedFile to write\n";

	my %resultHash;
	my $i = 0;

	my $npeak = 0;
	while (my $line = <$fin>)
	{
	    chomp $line;

	    print "$i ...\n" if $verbose && $i % 5000 == 0;
	    $i++;

	    my $peak = lineToBed ($line);
	    my $name = $peak->{'name'};
	    my $peakHeight = $peak->{'score'};

	    my ($peakId, $geneId) = split ("@@", $name);
		my ($leftBoundary, $rightBoundary, $leftHalfPH, $rightHalfPH);

		($peakId, $leftBoundary, $rightBoundary, $leftHalfPH, $rightHalfPH) = split ("#", $peakId) if $valleySeeking;

    	next unless exists $geneTagCountHash{$geneId};

    	my $expectedPeakHeight = $geneTagCountHash{$geneId}->{'PH0'}; #expected PH
    	next unless $expectedPeakHeight > 0; #otherwise no tag on the exonic region of the gene
    	next unless $peakHeight > $expectedPeakHeight;#otherwise it cannot be significant

    	my $geneSize = $geneTagCountHash{$geneId}->{'size'};

		#print "$i: $peakHeight, $expectedPeakHeight\n";
		if ($peakHeight >= 5000)
		{
			#print a warning
			print "$i: peak height > 5000: ", bedToLine ($peak), "\n";
		}


		my $pvalue;

		if ($maxPH > 0 && $peakHeight > $maxPH)
		{
			#peak height too large to evaluate scan statistics. flagged by nd in pvalue
			print "$i: peak height out of range: ", bedToLine ($peak), "\n" if $verbose;
			$pvalue = "nd";
			next if $skipOutofRangePeaks;
		}
		else
		{
    		$pvalue = exists $resultHash{$geneId}{$peakHeight} ? 
    			$resultHash{$geneId}{$peakHeight} : calcScanStatistic ($peakHeight, $expectedPeakHeight, $geneSize / $tagSize);

    		$resultHash{$geneId}{$peakHeight} = $pvalue unless exists $resultHash{$geneId}{$peakHeight};

    		$pvalue *= $effectiveGeneNum if $multiTestCorrection;
    		$pvalue = 1e-100 if $pvalue == 0;
			$pvalue = sprintf ("%.2e", $pvalue);
		}
		
    	$peak->{'name'} = $peakId . "[gene=$geneId][PH=$peakHeight][PH0=" . 
    		sprintf ("%.2f", $expectedPeakHeight) . "][P=$pvalue]";
	
		$peak->{'name'} .= "#" . join("#", $leftBoundary, $rightBoundary, $leftHalfPH, $rightHalfPH) if $valleySeeking;
    	$peak->{'score'} = $peakHeight; #-log ($pvalue) / log(10);
    
    	print $fout bedToLine ($peak), "\n" if $pvalue <= $pvalueThreshold;
		$npeak++;
	}

	close ($fin);
	close ($fout);
	Carp::croak "no significant peaks found\n" if $npeak == 0;
}


if ($maxGap >=0)
{
	print "merge peaks ...\n" if $verbose;
	my $tmpPeakMergedBedFile = "$cache/tmp.peak.sig.merge.bed";
	$cmd = "perl $cmdDir/bedUniq.pl $ssFlag -c maxScore -maxgap $maxGap $verboseFlag $tmpPeakBedFile $tmpPeakMergedBedFile";
	system ($cmd);
	$tmpPeakBedFile = $tmpPeakMergedBedFile;
}

print "generating output...\n" if $verbose;
#output, some formatting issues
my $fin;
my $fout;
my $foutPeak;
my $foutBoundary;
my $foutHalfPH;

open ($fin, "<$tmpPeakBedFile") || Carp::croak "cannot open file $tmpPeakBedFile to read\n";
open ($fout, ">$outBedFile") || Carp::croak "cannot open file $outBedFile to write\n";

if ($outBoundaryBedFile ne '')
{
	open ($foutBoundary, ">$outBoundaryBedFile") || Carp::croak "cannot open file $outBoundaryBedFile to write\n";
}

if ($outHalfPHBedFile ne '')
{
	open ($foutHalfPH, ">$outHalfPHBedFile") || Carp::croak "cannot open file $outHalfPHBedFile to write\n";
}

my $iter = 0;
while (my $line = <$fin>)
{
	chomp $line;
	next if $line=~/^\s*$/;

	$iter++;
	my $p = lineToBed ($line);
	my $peakId = $p->{'name'};
	my ($leftBoundary, $rightBoundary, $leftHalfPH, $rightHalfPH);
	($peakId, $leftBoundary, $rightBoundary, $leftHalfPH, $rightHalfPH) = split ("#", $peakId) if $valleySeeking;

	$p->{'name'} = sprintf("%s_%d", $prefix, $iter);
	if ($geneBedFile ne '')
	{
		$peakId =~/(\[\S*)$/;
		my $sig = $1;
		$p->{'name'} .= $sig;
	}
	print $fout bedToLine ($p), "\n";

	my %pCopy = %$p;
	if ($outBoundaryBedFile ne '')
	{
		$pCopy{'chromStart'} = $leftBoundary;
		$pCopy{'chromEnd'} = $rightBoundary - 1;
		print $foutBoundary bedToLine (\%pCopy), "\n";
	}

	if ($outHalfPHBedFile ne '')
	{
		$pCopy{'chromStart'} = $leftHalfPH;
		$pCopy{'chromEnd'} = $rightHalfPH - 1;
		print $foutHalfPH bedToLine (\%pCopy), "\n";
	}
}

close ($fin);
close ($fout);
close ($foutBoundary) if $outBoundaryBedFile ne '';
close ($foutHalfPH) if $outHalfPHBedFile ne '';

system ("rm -rf $cache") unless $keepCache;



sub extractPeaks
{
	my ($tags, $strand, $valleyDepth) = @_;
	
	my @tags;
	#extract the tags on the specified strand
	foreach my $t (@$tags)
    {
        $t->{"strand"} = '+' unless exists $t->{"strand"};
        next if $strand ne 'b' && $strand ne $t->{"strand"};
        push @tags, $t;
    }
	
	#sort the tags
	@tags = sort {$a->{'chromStart'} <=> $b->{'chromStart'}} @tags;	

	my @peaks;

	my @currBlock;
	my $currEnd = -1;
	foreach my $t (@tags)
	{
		if ($t->{'chromStart'} > $currEnd + 1)
		{
			if (@currBlock > 0)
			{
				my $ret = extractPeaksBlock (\@currBlock, $valleyDepth);
				push @peaks, @$ret;
			
				#start a new block
				@currBlock = ();
			}
		}
		push @currBlock, $t;
		$currEnd = $t->{'chromEnd'};
	}

	#the last block
	my $ret = extractPeaksBlock (\@currBlock, $valleyDepth);
	push @peaks, @$ret;

	return \@peaks;
}


#an upstream bounding local maixma is an upstream local maxima whose PH is greater than the current local maxima
#simiarly, a downstream bounding local maxima is a downstream local maxima whose PH is greater than the current local maxima

#a peak is a local maxima satisfying
#a local minima with PH< valley height can be found between the upstream bounding and the current local maxima
#a local minima with PH< valley height can be found between the downstream bounding and the current local maxima


sub extractPeaksBlock
{
	my ($tags, $valleyDepth) = @_;

	if (@$tags == 1)
	{
		my $t = $tags->[0];
		$t->{'leftHalfPH'} = $t->{'chromStart'};
		$t->{'rightHalfPH'} = $t->{'chromEnd'};
		$t->{'leftBoundary'} = $t->{'chromStart'};
		$t->{'rightBoundary'} = $t->{'chromEnd'};
		
		return $tags;
	}
	my $valleyHeight = 1 - $valleyDepth;


	#Carp::croak Dumper ($tags), "\n";

	#label local maxima and minima
	#assuming the tags are already sorted here
	my @maxima;
	my @minima;
	for (my $i = 0; $i < @$tags; $i++)
	{
		$tags->[$i]->{'idx'} = $i;
		if ($i == 0)
		{
			if ($tags->[$i]->{'score'} > $tags->[$i+1]->{'score'})
			{
				push @maxima, $i;
				$tags->[$i]->{'maxima'} = $#maxima;
			}
			elsif ($tags->[$i]->{'score'} < $tags->[$i+1]->{'score'})
			{
				push @minima, $i;
				$tags->[$i]->{'minima'} = $#minima;
			}
		}
		elsif ($i > 0 && $i < @$tags - 1)
		{
			if ($tags->[$i]->{'score'} > $tags->[$i+1]->{'score'} && $tags->[$i]->{'score'} >= $tags->[$i-1]->{'score'})
			{
				push @maxima, $i;
				$tags->[$i]->{'maxima'} = $#maxima;
			}
			elsif ($tags->[$i]->{'score'} < $tags->[$i+1]->{'score'} && $tags->[$i]->{'score'} <= $tags->[$i-1]->{'score'})
			{
				push @minima, $i;
				$tags->[$i]->{'minima'} = $#minima;
			}
		}
		else
		{
			if ($tags->[$i]->{'score'} >= $tags->[$i-1]->{'score'})
			{
				push @maxima, $i;
				$tags->[$i]->{'maxima'} = $#maxima;
			}
			elsif ($tags->[$i]->{'score'} <= $tags->[$i-1]->{'score'})
			{
				push @minima, $i;
				$tags->[$i]->{'minima'} = $#minima;
			}
		}
	}
	#Carp::croak Dumper ($tags), "\n";
	
	my @peaks;
	#find the closest neighbor maxima which has bigger peak height
	for (my $i = 0; $i < @maxima; $i++)
	{
		#this is naive implementation, could improve speed further
		my $t1 = $tags->[$maxima[$i]];
		for (my $j = $i - 1; $j>= 0; $j--)
		{
			my $t2 = $tags->[$maxima[$j]];

			#to break ties
			#if both PH and size are equal, we go further on the right, but not on the left
			if ($t1->{'score'} < $t2->{'score'} || ($t1->{'score'} == $t2->{'score'} && $t1->{'chromEnd'} - $t1->{'chromStart'} < $t2->{'chromEnd'} - $t2->{'chromStart'}))
			{
				$t1->{'prev'} = $maxima[$j];
				last;
			}
		}

		for (my $j = $i + 1; $j < @maxima; $j++)
		{
			my $t2 = $tags->[$maxima[$j]];
			#if both PH and size are equal, we go further on the right, but not on the left
			if ($t1->{'score'} < $t2->{'score'} || ($t1->{'score'} == $t2->{'score'} && $t1->{'chromEnd'} - $t1->{'chromStart'} <= $t2->{'chromEnd'} - $t2->{'chromStart'}))
            {
                $t1->{'next'} = $maxima[$j];
                last;
            }
		}
	}
	#Carp::croak Dumper ($tags), "\n";

	#find if a valley exists surrounding each local maxima
	for (my $i = 0; $i < @maxima; $i++)
    {
		my $p1 = $tags->[$maxima[$i]];
		
		my $leftValleyFound = 0;
		my $rightValleyFound = 0;

		if (exists $p1->{'prev'})
		{
			my $prev = $p1->{'prev'};
			for (my $j = $maxima[$i] - 1; $j >= $prev+1; $j--)
			{
				my $v = $tags->[$j];
				next unless exists $v->{'minima'};
				if ($v->{'score'} / $p1->{'score'} <= $valleyHeight)
				{
					$leftValleyFound = 1;
					last;
				}
			}
		}
		else
		{
			$leftValleyFound = 1;
		}

		if (exists $p1->{'next'})
		{
			my $next = $p1->{'next'};
			for (my $j = $maxima[$i] + 1; $j < $next; $j++)
			{
				my $v = $tags->[$j];
				next unless exists $v->{'minima'};
				if ($v->{'score'} / $p1->{'score'} <= $valleyHeight)
                {
                    $rightValleyFound = 1;
                    last;
                }
			}
		}
		else
		{
			$rightValleyFound = 1;
		}

		push @peaks, $p1 if $leftValleyFound && $rightValleyFound;
	}
	
	for (my $i = 0; $i < @peaks; $i++)
	{
		my $p1 = $peaks[$i];
		#search peak boundary and half PH boundary
		my $leftBoundaryTag = $p1;
		my $leftHalfPHTag = 0;
			
		my $left = exists $p1->{'prev'} ? $p1->{'prev'} + 1 : 0;
		if ($i > 0)
		{
			$left = $peaks[$i-1]->{'idx'} + 1 if $left < $peaks[$i-1]->{'idx'} + 1;
		}
		for (my $j = $p1->{'idx'} - 1; $j >= $left; $j--)
		{
			my $v = $tags->[$j];
			if ($v->{'score'} / $p1->{'score'} <= 0.5)
			{
				#the first time the coverage reduced to 1/2PH
				$leftHalfPHTag = $v unless $leftHalfPHTag;
			}
				
			if ($leftBoundaryTag->{'score'} > $v->{'score'})
			{
				$leftBoundaryTag = $v;
			}
		}
		$leftHalfPHTag = $leftBoundaryTag unless $leftHalfPHTag;

		my $rightBoundaryTag = $p1;
		my $rightHalfPHTag = 0;
		my $right = exists $p1->{'next'} ? $p1->{'next'} - 1 : $#$tags;
		if ($i < @peaks - 1)
		{
			$right = $peaks[$i+1]->{'idx'} - 1 if $right > $peaks[$i+1]->{'idx'} - 1;
		}

		for (my $j = $p1->{'idx'} + 1; $j <= $right; $j++)
		{
			my $v = $tags->[$j];
			if ($v->{'score'} / $p1->{'score'} <= 0.5)
            {
				$rightHalfPHTag = $v unless $rightHalfPHTag;
			}
			
			if ($rightBoundaryTag->{'score'} > $v->{'score'})
			{
				$rightBoundaryTag = $v;
			}
		}
		$rightHalfPHTag = $rightBoundaryTag unless $rightHalfPHTag;
		$p1->{'leftHalfPH'} = $leftHalfPHTag->{'idx'} < $p1->{'idx'} ? ($leftHalfPHTag->{'chromEnd'} + 1) : $leftHalfPHTag->{'chromStart'};
		$p1->{'rightHalfPH'} = $rightHalfPHTag->{'idx'} > $p1->{'idx'} ? ($rightHalfPHTag->{'chromStart'} - 1) : $rightHalfPHTag->{'chromEnd'};
		$p1->{'leftBoundary'} = $leftBoundaryTag->{'idx'} < $p1->{'idx'} ? $leftBoundaryTag->{'chromEnd'} + 1 : $leftBoundaryTag->{'chromStart'};
		$p1->{'rightBoundary'} = $rightBoundaryTag->{'idx'} > $p1->{'idx'} ? $rightBoundaryTag->{'chromStart'} - 1 : $rightBoundaryTag->{'chromEnd'};
		
	}

=debug
	if ($tags->[0]->{'chrom'} eq 'chr1' && $tags->[0]->{'strand'} eq '+' && $tags->[0]->{'chromStart'} > 66439247 && $tags->[$#$tags]->{'chromEnd'} < 66439590)
	{
		print Dumper ($tags), "\n";
		print "maxima=", join ("\t", @maxima), "\n";
		print "minima=", join ("\t", @minima), "\n";
		Carp::croak Dumper (\@peaks), "\n";
	}
=cut
	return \@peaks;
}

sub getLocationInfo
{
    my ($conf, $dbkey, $analysis) = @_;

    my $fin;

    open ($fin, "<$conf") || Carp::croak "cannot open file $conf\n";

    my %ret;

    while (my $line = <$fin>)
    {
        chomp $line;
        next if $line =~/^\s*$/;
        next if $line =~/^\#/;

        my ($db, $ana, $path, $type) = split (/\s+/, $line);

        $path = "$cmdDir/$path" unless $path=~/^\//;
        #if a relative path is provided, we assume the annotation file is located in the same folder as the script
        #fixed by CZ, 07/31/2016

        $type = $ana unless $type;

        if ($db eq $dbkey && $ana eq $analysis)
        {
            Carp::croak "$path does not exist\n" unless -f $path;
            $ret{$type} = $path;
        }
    }
    close ($fin);
    return \%ret;
}

