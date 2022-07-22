#!/usr/bin/env perl

#this algorithm is efficient when the reads are the same (or similar ) sizes
#it calculate tag numbers either in fixed size moving windows or in specified regions
#for the latter, i assume the latter is not too large, so that it can be stored in the memory


use strict;
use warnings;

use Getopt::Long;
use Carp;
use File::Basename;
use Data::Dumper;

use MyConfig;
use Common;
use Bed;

STDOUT->flush();
STDERR->flush();


my $prog = basename ($0);

my $bigFile = 0;
my $minBlockSize = 2000000;
my $ext5 = 0;
my $ext3 = 0;
my $weight = 0;
my $weightAvg = 0;
my $separateStrand = 0;
my $chromLenFile = "";

my $regionFile = "";
my $exact = 0;
my $windowSize = 100;
my $stepSize = 20;
my $normalizationMethod = 'none'; # or rpkm, multiply={1.3}
my $normalizationFactor = 1;



my $verbose = 0;
my $outputFormat = 'bed'; #wig or bedgraph or sgr
my $noZero = 0;
my $trackName = '';

my $cache = getDefaultCache ($prog);
my $keepCache = 0;

my @ARGV0 = @ARGV;

GetOptions (
		'big'=>\$bigFile,
		'minBlockSize:i'=>\$minBlockSize,
		'chrLen:s'=>\$chromLenFile,
		'ext5:i'=>\$ext5,
		'ext3:i'=>\$ext3,
		'weight'=>\$weight,
		'weight-avg'=>\$weightAvg,
		'region:s'=>\$regionFile,
		'exact'=>\$exact,
		'normalize:s'=>\$normalizationMethod,
		'w|window-size:i'=>\$windowSize,
		's|step-size:i'=>\$stepSize,
		'ss|sep-strand'=>\$separateStrand,
		'n|name:s'=>\$trackName,
		'of:s'=>\$outputFormat,
		'nz'=>\$noZero,
		'c|cache:s'=>\$cache,
		'keep-cache'=>\$keepCache,
		'v|verbose'=>\$verbose);


if (@ARGV != 2 && @ARGV != 3)
{
	print "Count tag numbers in sliding window or particular regions or each nucleotide position\n";
	print "Usage: $prog [options] <in.bed> <profile.out> [profile.out.2]\n";
	print " <in.bed> : gzip compressed Bed file with .gz extension is allowed\n";
	print "          : use \"-\" for stdin\n"; 
	print " <profile.out> [profile.out.2] : for wig format output, specify two files to separate the two strands\n";
	print "OPTIONS:\n";
	
	print "\n[Input options]\n";
	print " -big                : set when the file is big\n";
	print " -minBlockSize  [int]: minimim number of lines to read in each block for a big file ($minBlockSize)\n";
	print " -weight             : weight counts according to the score of each tag\n";
	print " -weight-avg         : weight average the score of each tag\n";
	print " -ss                 : separate strand (off)\n";
	print " -ext5          [int]: extension of tags at the 5' end\n";
	print " -ext3          [int]: extension of tags at the 3' end\n";
	print " -chrLen     [string]: chrom len file\n";
	
	print "\n[Profile options]\n";
	print " -region       [file]: a bed file with regions to count tag numbers. if not specified, count in moving windows\n";
	print " -exact              : exact count at each nucleotide\n";
	print " -w             [int]: window size ($windowSize)\n";
	print " -s             [int]: step size ($stepSize)\n";
	print " --normalize [string]: normalization ([none]|rpkm|multiply={1.3})\n";
	
	#print " --normalizing-factor: normalizing factor ($normalizingFactor)\n";
	print "\n[Output options]\n";
	print " -of         [string]: output format ([bed]|bedgraph|sgr)\n";
	print " -nz                 : don't print zeroes (work for sgr and bed)\n";
	print " -n          [string]: track name ($trackName)\n";
	print " -c          [string]: cache dir ($cache)\n";
	print " --keep-cache        : keep cache when the job done\n";
	print " -v                  : verbose\n";
	exit (1);
}

print "CMD=$prog ", join(" ", @ARGV0), "\n" if $verbose;

my ($tagBedFile, $outProfileFile, $outProfileFile2) = @ARGV;

$outProfileFile2 = "" unless $outProfileFile2;

$outputFormat = "bedgraph" if $outputFormat eq 'wig';

#check parameters
if ($outProfileFile2)
{
	Carp::croak "output format must be bedgraph and separate strand must be specified\n"
	unless $outputFormat eq 'bedgraph' && $separateStrand == 1;
}


if ($outputFormat eq 'bedgraph' && $separateStrand == 1 && $outProfileFile2 eq '')
{
	#if we are going to put two tracks in the same file, we have to have trackName
	$trackName = "profile" unless $trackName;
}


my $profileType  = "window";

if ($exact)
{
	$profileType = 'exact';
}
elsif (-f $regionFile)
{
	$profileType = 'region';
}

if ($weight && $weightAvg)
{
	Carp::croak "-weight and -weight-avg cannot be set at the same time\n";
}

if ($regionFile ne '' && $exact)
{
	Carp::croak "--exact and -region cannot be set at the same time\n";
}

if ($normalizationMethod=~/^multiply/)
{
	$normalizationMethod=~/^multiply=(\S*)$/;
	$normalizationFactor = $1;

	Carp::croak "normalization factor must be positive\n" unless $normalizationFactor > 0;
}
elsif ($normalizationMethod eq 'rpkm')
{
	Carp::croak "-weight-avg does not work for rpkm\n" if $weightAvg;
}
elsif ($normalizationMethod ne '' && $normalizationMethod ne 'none' && $normalizationMethod ne 'rpkm')
{
	Carp::croak "incorrect option --normalize $normalizationMethod\n";
}

if ($profileType eq 'window')
{
	Carp::croak "no chrLen.txt file specified\n" unless -f $chromLenFile;
}


die "$cache already exist\n" if (-d $cache);
die "$cache already exist\n" if (-f $cache);

system ("mkdir $cache");

my $currOutputLineIter = 0; #to keep track how many lines of output, used in the name column of bed files in some cases


#loading input data
my $fin;
my %chromLen;
if (-f $chromLenFile)
{
	print "reading chrom length from $chromLenFile ...\n" if $verbose;
	open ($fin, "<$chromLenFile") || Carp::croak "can not open file $chromLenFile to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		my ($chrom, $len) = split (/\t/, $line);
		$chromLen{$chrom} = $len;
	}

	close ($fin);
}


print "reading tags from $tagBedFile ...\n" if $verbose;

my %tagCount;
if ($bigFile)
{
	my $ret = splitBedFileByChrom ($tagBedFile, $cache, v=>$verbose, "sort"=>1);
	%tagCount = %$ret;
}
else
{
	my $tags = readBedFile ($tagBedFile, $verbose);
	foreach my $t (@$tags)
	{
		my $chrom = $t->{"chrom"};
		push @{$tagCount{$chrom}}, $t;
	}		
}

print "get tag count broken down into chromosomes ...\n" if $verbose;

foreach my $chrom (sort keys %tagCount)
{
	
	my $n = $tagCount{$chrom};
	$n = ref($n) eq 'HASH' ? $n = $n->{'n'} : @$n;
	print "$chrom : $n tags\n" if $verbose;
}

#read regions to profile if necessary
my %regionHash= ();
my $regions = [];
if ($profileType eq 'region')
{

	print "reading regions to count tag numbers from $regionFile ...\n" if $verbose;
	
	$regions = readBedFile ($regionFile);
	my $n = @$regions;
	print "$n regions loaded\n" if $verbose;

	foreach my $r (@$regions)
	{
		my $chrom = $r->{"chrom"};
		$r->{"score"} = 0; #initialize the score

		$r->{"name"} = $r->{"chrom"} . ":" . ($r->{"chromStart"} + 1) . "-" . ($r->{"chromEnd"} + 1) unless exists $r->{"name"};

		if ($separateStrand)
		{
			Carp::croak "no strand information in ", Dumper ($r), "\n" unless exists $r->{"strand"};
			my $strand = $r->{"strand"};
			
			push @{$regionHash{$strand}->{$chrom}}, $r;
		}
		else
		{
			push @{$regionHash{'b'}->{$chrom}}, $r;
		}
	}
}

#calculate normalization factor if necessary
if ($normalizationMethod eq 'rpkm')
{
	print "count the total number of tags and average tag size ...\n" if $verbose;
	my $totalTagNum = 0;
	my $totalNucleotideCount = 0;

	foreach my $chrom (sort keys %tagCount)
	{	
		print "processing tags on chrom $chrom ...\n" if $verbose;

		if ($bigFile)
		{
			my $tmpFile = $tagCount{$chrom}->{'f'};
			my $fin;
			open ($fin, "<$tmpFile") || Carp::croak "cannot open file $tmpFile to read\n";
			my $iter = 0;
			while (my $t = readNextBedLine ($fin))
			{
				print "$iter ...\n" if $verbose && $iter % 100000 == 0;
				$iter++;
				my $tagLen = exists $t->{'blockSizes'} ? sum ($t->{'blockSizes'}) : ($t->{'chromEnd'} - $t->{'chromStart'} + 1);
				my $copyNum = $weight ? $t->{'score'} : 1;
				$totalTagNum += $copyNum;
				$totalNucleotideCount += $tagLen * $copyNum;
			}
			close ($fin);
		}
		else
		{
			my $tags = $tagCount{$chrom};
			foreach my $t (@$tags)
			{
				my $tagLen = exists $t->{'blockSizes'} ? sum ($t->{'blockSizes'}) : ($t->{'chromEnd'} - $t->{'chromStart'} + 1);
				my $copyNum = $weight ? $t->{'score'} : 1;

				$totalTagNum += $copyNum;
				$totalNucleotideCount += $tagLen * $copyNum;
			}
		}
	}

	print "$totalNucleotideCount nucleotides in $totalTagNum tags found\n";

	#$averageTagLen = $totalNucleotideCount / $totalTagNum;
	if ($profileType eq 'exact')
	{
		$normalizationFactor = 1e9 / $totalNucleotideCount;
	}
	elsif ($profileType eq 'region')
	{
		$normalizationFactor = 1e9 / $totalTagNum; #note that we cannot divide the number by the size of each region now
	}
	else
	{
		$normalizationFactor = 1e9 / $totalTagNum / $windowSize;
	}

	print "updated normalization factor: $normalizationFactor\n" if $verbose;
}


my @strand = ('b');
@strand = qw(+ -) if $separateStrand;

my ($fout, $fout2);
open ($fout, ">$outProfileFile") || Carp::croak "can not open file $outProfileFile to write\n";

if ($outProfileFile2)
{
	open ($fout2, ">$outProfileFile2") || Carp::croak "can not open file $outProfileFile2 to write\n";
}


foreach my $s (@strand)
{
	
	print "#processing strand $s ...\n" if $verbose;
	if ($profileType eq 'region')
	{
		next unless exists $regionHash{$s};
	}

	if ($outputFormat eq 'bedgraph' && $trackName ne '' && (not (-f $regionFile)))
	{
		#we need track headers
		if ($s eq 'b')
		{
			print $fout "\n\ntrack type=bedGraph name=\"$trackName\" autoScale=on maxHeightPixels=128:36:16\n"; # windowingFunction=mean\n";
		}
		else
		{
			if ($s eq '-' && $outProfileFile2 ne '')
			{
				print $fout2 "\n\ntrack type=bedGraph name=\"$trackName($s)\" autoScale=on maxHeightPixels=128:36:16\n";# windowingFunction=mean\n";
			}
			else
			{
				print $fout "\n\ntrack type=bedGraph name=\"$trackName($s)\" autoScale=on maxHeightPixels=128:36:16\n";# windowingFunction=mean\n";
			}
		}
	}
	
	my @chroms = $profileType eq 'window' ? sort keys %chromLen : sort keys %tagCount;
	#for the sliding window mode, we generate output even when there are no tags
	
	foreach my $chrom (@chroms)
	{
		if ($profileType eq 'region')
		{
			next unless exists $regionHash{$s}->{$chrom};
		}

		print "processing chromsome $chrom ...\n" if $verbose;
		
		$currOutputLineIter = 0; #reset
		if (exists $tagCount{$chrom})
		{

			if ($bigFile)
			{

				my $tmpFile = $tagCount{$chrom}->{'f'};
				my $fin;
				open ($fin, "<$tmpFile") || Carp::croak "cannot open file $tmpFile to read\n";
				
				my $iter = 0;
				my @profile;
				
				while (my $block = readNextBedBlock ($fin, max (1, $ext5 + $ext3 + 1), 0, minBlockSize=> $minBlockSize))
				{
					my $n = @$block;
					print "block $iter : $n tags\n" if $verbose;
					$iter++;

					if ($profileType eq 'window')
					{
						#window profile, do not print here, but sum them up
						my $ret = getProfile ($block, $chrom, $s, $fout, $fout2, "print"=> 0);
						#Carp::croak Dumper ($ret), "\n";
						my $n = @$ret;
						map
						{
							$profile[$_]->{'n'} += $ret->[$_]->{'n'}; 
							$profile[$_]->{'n2'} += $ret->[$_]->{'n2'};
						} (0..($n-1));
					}
					else
					{
						#exact profile will print output immediately, but region profile will print at the end
						getProfile ($block, $chrom, $s, $fout, $fout2);
					}
					
				}
				close ($fin);

				#output window profile
				printWindowProfile ($chrom, \@profile, $s, $fout, $fout2) if $profileType eq 'window';
			}
			else
			{
				my $tags = $tagCount{$chrom};
				my $n = @$tags;
				print "$n tags loaded on chromosome $chrom\n" if $verbose;
				#print everything immediately
				getProfile ($tags, $chrom, $s, $fout, $fout2);
			}
		}
		else
		{
			my $tags = [];
			#only for windows profile, the other types of profile will do nothing
			getProfile ($tags, $chrom, $s, $fout, $fout2);
		}
	}#chrom
}#strand


close ($fout);
close ($fout2) if ($outProfileFile2 ne '');


if ($profileType eq 'region')
{
	print "writing profile to $outProfileFile ...\n" if $verbose;
	
	foreach my $r (@$regions)
	{
		if (exists $r->{'n'})
		{
			my $score = ($weightAvg && $r->{'n2'}) ? $r->{'n'} / $r->{'n2'} : $r->{'n'}; 
			$score *= $normalizationFactor;
			$score /= ($r->{'chromEnd'} - $r->{'chromStart'} + 1) if $normalizationMethod eq 'rpkm';
			$r->{'score'} = $score;
		}
	}

	if ($outputFormat eq 'bedgraph')
	{
		my $fout;
		open ($fout, ">$outProfileFile") || Carp::croak "cannot openfile $outProfileFile to write\n";
		print $fout "\n\ntrack type=bedGraph name=\"$trackName\" autoScale=on maxHeightPixels=128:36:16\n" if $trackName ne '';
		foreach my $r (@$regions)
		{
			print $fout join ("\t", $r->{'chrom'}, $r->{'chromStart'}, $r->{'chromEnd'} + 1, $r->{'score'}), "\n";
		}
		close ($fout);
	}
	else
	{
		writeBedFile ($regions, $outProfileFile);
	}
}
system ("rm -rf $cache") unless $keepCache;


#this subroutine used global variable, so should never be moved to the library
sub getProfile
{
	my ($tags, $chrom, $strand, $fout, $fout2, %params) = @_;

	my $print = exists $params{'print'} ? $params{'print'} : 1; #only for window profile
	
	
	#get tags on specific strand
	my @tags;
	foreach my $t (@$tags)
	{
		Carp:croak "cannot extend reads without strand information\n" if ($ext5 != 0 || $ext3 != 0) & not exists $t->{"strand"};

		$t->{"strand"} = '+' unless exists $t->{"strand"};
			
		next if $strand ne 'b' && $strand ne $t->{"strand"};
			
		#extend tags if necessary
		if ($t->{"strand"} eq '+')
		{
			$t->{"chromStart"} -= $ext5;
			$t->{"chromEnd"} += $ext3;
		}
		else
		{
			$t->{"chromStart"} -= $ext3;
			$t->{"chromEnd"} += $ext5;
		}
	
	
		push @tags, $t;
	}

	#my $chrom = $tags->[0]->{'chrom'};

	if ($profileType eq 'exact')
	{
		my $profile = getExactProfile (\@tags);
		#print Dumper ($profile), "\n";	
		for (my $i = 0; $i < @$profile - 1; $i++)
		{
			my $chromStart = $profile->[$i]->[0];
			my $chromEnd = $profile->[$i+1]->[0] -1;
			my $score = $profile->[$i]->[1];

			$score *= $normalizationFactor;
			$score = sprintf ("%.2f", $score) if $score  * 100 - int ($score * 100) != 0;


			if ($outputFormat eq 'bedgraph')
			{
				if ($outProfileFile2 ne '' && $strand eq '-')
				{
					#we alerady check the consistency of parameters, of the existence of the second ouput file should be surfice to print into a separate file handle
					print $fout2 join ("\t", $chrom, $chromStart, $chromEnd+1, $score), "\n" if $score > 0;
				}
				else
				{
					print $fout join ("\t", $chrom, $chromStart, $chromEnd+1, $score), "\n" if $score > 0;
				}
			}
			elsif ($outputFormat eq 'bed')
			{
				if ($strand eq 'b')
				{
					print $fout join ("\t", $chrom, $chromStart, $chromEnd+1, $i + $currOutputLineIter, $score), "\n" if $score > 0;
				}
				else
				{
					print $fout join ("\t", $chrom, $chromStart, $chromEnd+1, $i + $currOutputLineIter, $score, $strand), "\n" if $score > 0;
				}
			}
			else
			{
				Carp:croak "output format $outputFormat is not implemented for the exact mode\n";
			}
		}
		$currOutputLineIter += @$profile;
		return $profile;
	}
	elsif ($profileType eq 'region')
	{
		return getRegionProfile (\@tags, $regionHash{$strand}->{$chrom});
=annot		
		for (my $i = 0; $i < @$profile; $i++)
		{
			#next if $profile->[$i]->{"n"} == 0 && $noZero;
			my $r = $profile->[$i]->{"r"};
			$r->{"score"} = $profile->[$i]->{"n"} * $normalizationFactor;
			$r->{'score'} = $r->{'score'} / ($r->{'chromEnd'} - $r->{'chromStart'} + 1) if $normalizationMethod eq 'rpkm';
		}
=cut
			#we do not output regions now
	}
	else
	{
		#window profile
		my $profile = getWindowProfile (\@tags, $windowSize, $stepSize, $chromLen{$chrom});

		#Carp::croak "profile=", Dumper ($profile), "\n"; # if ($chrom eq 'chr13_random');
		printWindowProfile ($chrom, $profile, $strand, $fout, $fout2) if $print;
		
		return $profile;
	}
}

#this subroutine used global variable, so should never be moved to the library
sub printWindowProfile
{
	my ($chrom, $profile, $strand, $fout, $fout2) = @_;

	my @p;
	for (my $i = 0; $i < @$profile; $i++)
	{
		my $score = $profile->[$i]->{'n'};
		$score /= $profile->[$i]->{'n2'} if ($weightAvg && $profile->[$i]->{'n2'});
		push @p, $score * $normalizationFactor;
	}


	if ($outputFormat eq 'sgr')
	{
		for (my $i = 0; $i < @$profile; $i++)
		{
			next if $p[$i] == 0 && $noZero;

			my $chromStart = $i * $stepSize;
			my $pos = $chromStart + int ($windowSize / 2);

			if ($strand eq 'b')
			{
				print $fout join ("\t", $chrom, $pos, $p[$i]), "\n";
			}
			else
			{
				print $fout join ("\t", $chrom, $pos, $p[$i], $strand), "\n";
			}
		}
	}
	elsif ($outputFormat eq 'bedgraph')# wig
	{
		my $start = int ($windowSize / 2);
		if ($strand eq '-' && $outProfileFile2 ne '')
		{
			print $fout2 "fixedStep chrom=$chrom start=$start step=$stepSize\n";
			print $fout2 join ("\n", @p), "\n";
		}
		else
		{
			print $fout "fixedStep chrom=$chrom start=$start step=$stepSize\n";
			print $fout join ("\n", @p), "\n";
		}
	}
	else # bed
	{
		for (my $i = 0; $i < @$profile; $i++)
		{
			next if $p[$i] == 0 && $noZero;
			my $chromStart = $i * $stepSize;
			my $chromEnd = $chromStart + $windowSize;
			if ($strand eq 'b')
			{
				print $fout join ("\t", $chrom, $chromStart, $chromEnd, "w" . ($i+$currOutputLineIter), $p[$i]), "\n";
			}
			else
			{
				print $fout join ("\t", $chrom, $chromStart, $chromEnd, "w" . ($i+$currOutputLineIter), $p[$i], $strand), "\n";
			}
		}
		$currOutputLineIter += @$profile;
	}
}

#this subroutine counts the number of tags in the window directly
#no update
sub getWindowProfile
{
	#tags are assumed to be on the proper strand of the same chromosome, but unsorted
	#
	my ($tags, $windowSize, $stepSize, $chrLen) = @_;
	
	#sort all tags according to the start coordinates
	#
	my @tags = sort {$a->{"chromStart"} <=> $b->{"chromStart"}} @$tags;
	
	#we also need to get a bound of the end coordinates
	#chromEndExt is >= chromEnd & chromEndExt [i-1] <= chromEndExt [i]
	
	my $chromEndExt = 0;
	foreach my $t (@tags)
	{
		$chromEndExt = $t->{"chromEnd"} if $t->{"chromEnd"} > $chromEndExt;
		$t->{"chromEndExt"} = $chromEndExt;
	}

	$tags = \@tags;
	
	my $firstTagIdx = 0; #the first tag that overlap with or on the right of the current window

	my @profile;
	for (my $windowStart = 0; $windowStart + $windowSize  <= $chrLen; $windowStart += $stepSize)
	{
		my $windowEnd = $windowStart + $windowSize - 1;
		my $tagNum = 0;
		my $tagNumNoWeight = 0;
		
		while ($firstTagIdx < @$tags && $tags->[$firstTagIdx]->{"chromEndExt"} < $windowStart)
		{
			$firstTagIdx++;
		}
		
		if ($firstTagIdx >= @$tags)
		{
			push @profile, {n=>0, n2=>0};
			next;
		}
		
		my $i = $firstTagIdx;
		
		while ($i < @$tags && $tags->[$i]->{"chromStart"} <= $windowEnd)
		{
			if ($tags->[$i]->{"chromEnd"} >= $windowStart)
			{
				my $overlap = 1;
				if (exists $tags->[$i]->{"blockCount"})
				{
					#check if the region is in the intron of the tags
					for (my $k = 1; $k < $tags->[$i]->{"blockCount"}; $k++)
					{
						my $intronStart = $tags->[$i]->{"chromStart"} + $tags->[$i]->{"blockStarts"}->[$k-1] + $tags->[$i]->{"blockSizes"}->[$k-1];
						my $intronEnd = $tags->[$i]->{"chromStart"} + $tags->[$i]->{"blockStarts"}->[$k] - 1;

						if ($intronStart <= $windowStart && $intronEnd >= $windowEnd)
						{
							$overlap = 0;
							last;
						}
					}
				}
				if ($overlap)
				{
					if ($weight || $weightAvg)
					{
						$tagNum += $tags->[$i]->{"score"};
					}
					else
					{
						$tagNum++; # if $tags->[$i]->{"chromEnd"} >= $windowStart;
					}
					$tagNumNoWeight++;
				}
			}

			$i++;
		}

		push @profile, {n=>$tagNum, n2=>$tagNumNoWeight};
=annot
		if ($weightAvg && $tagNumNoWeight)
		{
			push @profile, $tagNum / $tagNumNoWeight;
		}
		else
		{
			push @profile, $tagNum;
		}
=cut
	}

	return \@profile;
}


sub getRegionProfile
{
	my ($tags, $regions) = @_;
	#all tags and regions are assumed to be on the proper strand of the same chromosome
	#so strand information was not explicitly considered in the function
	#


	#sort all tags according to the start coordinates
	#
	my @tags = sort {$a->{"chromStart"} <=> $b->{"chromStart"}} @$tags;
	
	#we also need to get a bound of the end coordinates
	#chromEndExt is >= chromEnd & chromEndExt [i-1] <= chromEndExt [i]
	
	my $chromEndExt = 0;
	foreach my $t (@tags)
	{
		$chromEndExt = $t->{"chromEnd"} if $t->{"chromEnd"} > $chromEndExt;
		$t->{"chromEndExt"} = $chromEndExt;
	}

	$tags = \@tags;

	#sort regions
	my @regions = sort {$a->{"chromStart"} <=> $b->{"chromStart"}} @$regions;
	$regions = \@regions;

	
	my $firstTagIdx = 0; #the first tag that overlap with or on the right of the current window
	my @profile;
	foreach my $r (@$regions)
	{
		my $windowStart = $r->{"chromStart"};
		my $windowEnd = $r->{"chromEnd"};
		
		my $tagNum = 0;
		my $tagNumNoWeight = 0;

		while ($firstTagIdx < @$tags && $tags->[$firstTagIdx]->{"chromEndExt"} < $windowStart)
		{
			$firstTagIdx++;
		}
		
		my $i = $firstTagIdx;
		
		while ($i < @$tags && $tags->[$i]->{"chromStart"} <= $windowEnd)
		{
			if ($tags->[$i]->{"chromEnd"} >= $windowStart)
			{
				my $overlap = 1;
				if (exists $tags->[$i]->{"blockCount"})
				{
					#check if the region is in the intron of the tags
					for (my $k = 1; $k < $tags->[$i]->{"blockCount"}; $k++)
					{
						my $intronStart = $tags->[$i]->{"chromStart"} + $tags->[$i]->{"blockStarts"}->[$k-1] + $tags->[$i]->{"blockSizes"}->[$k-1];
						my $intronEnd = $tags->[$i]->{"chromStart"} + $tags->[$i]->{"blockStarts"}->[$k] - 1;

						if ($intronStart <= $windowStart && $intronEnd >= $windowEnd)
						{
							$overlap = 0;
							last;
						}
					}
				}

				if ($overlap)
				{
	
					if ($weight || $weightAvg)
					{
						$tagNum += $tags->[$i]->{"score"};
					}
					else
					{
						$tagNum++; # if $tags->[$i]->{"chromEnd"} >= $windowStart;
					}
					$tagNumNoWeight++;
				}
			}
			$i++;
		}

		$r->{'n'} += $tagNum;
		$r->{'n2'} += $tagNumNoWeight;
=annot

		if ($weightAvg && $tagNumNoWeight)
		{
			push @profile, {r=>$r, n=>$tagNum / $tagNumNoWeight};
		}
		else
		{
			push @profile, {r=>$r, n=>$tagNum};
		}
=cut
	}
	#return \@profile;
}


sub getExactProfile
{
	my ($tags) = @_;
	
	#tags are on the proper strand of the same chromosome
	
	my %tagCumulativeHash;
	my %tagCumulativeHashNoWeight;

	for my $t (@$tags)
	{
		my @blockStarts;
		my @blockSizes;
		
		if (exists $t->{"blockCount"})
		{
			@blockStarts = @{$t->{"blockStarts"}};
			@blockSizes = @{$t->{"blockSizes"}};
		}
		else
		{
			push @blockStarts, 0;
			push @blockSizes, $t->{"chromEnd"} - $t->{"chromStart"} + 1;
		}
		
		for (my $i = 0; $i < @blockSizes; $i++)
		{
			my $chromStart = $t->{"chromStart"} + $blockStarts[$i];
			my $chromEnd = $chromStart + $blockSizes[$i] - 1; #$t->{"chromEnd"};

			$tagCumulativeHash{$chromStart} = 0 unless exists $tagCumulativeHash{$chromStart};
			$tagCumulativeHash{$chromEnd+1} = 0 unless exists $tagCumulativeHash{$chromEnd+1};
		
			$tagCumulativeHashNoWeight{$chromStart} = 0 unless exists $tagCumulativeHashNoWeight{$chromStart};
			$tagCumulativeHashNoWeight{$chromEnd+1} = 0 unless exists $tagCumulativeHashNoWeight{$chromEnd+1};
		
			if ($weight || $weightAvg)
			{
				$tagCumulativeHash{$chromStart} += $t->{"score"};
				$tagCumulativeHash{$chromEnd+1} -= $t->{"score"};
			}
			else
			{
				$tagCumulativeHash{$chromStart}++;
				$tagCumulativeHash{$chromEnd+1}--;
			}

			$tagCumulativeHashNoWeight{$chromStart}++;
			$tagCumulativeHashNoWeight{$chromEnd+1}--;
		}
	}
	
	my @profile;

	my $count = 0;
	my $countNoWeight = 0;

	foreach my $pos (sort {$a <=> $b} keys %tagCumulativeHash)
	{
		$count += $tagCumulativeHash{$pos};
		$countNoWeight += $tagCumulativeHashNoWeight{$pos};

		Carp::croak "tag count is less than zero at $pos\n" if $countNoWeight < 0;

		my $intensity = $weightAvg && $countNoWeight ? ($count / $countNoWeight) : $count;
	
		if (@profile == 0 || Common::ABS($profile[$#profile]->[1] - $intensity) > 1e-6)
		{
			#to be more compact, we do not add new point if the intensity does not change !!
			push @profile, [$pos, $intensity];
		}

		#if ($weightAvg && $countNoWeight)
		#{
		#	push @profile, [$pos, $count / $countNoWeight];
		#}
		#else
		#{
		#	push @profile, [$pos, $count];
		#}
	}

	#Carp::croak Dumper (\@profile), "\n";

	if (@profile > 0)
	{
		my $lastNode = $#profile;
		Carp::croak "tag count is not zero at last node:", Dumper (\@profile), "\n" if Common::ABS ($profile[$lastNode][1]) > 1e-4;
	}
	return \@profile;
}

