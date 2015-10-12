#!/usr/bin/perl -w

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Carp;
use Data::Dumper;

use Bed;
use Common;


my $verbose = 0;
my $prog = basename ($0);
my $sameStrand = 0;
#my $step = 1;
my $noMatchScore = -1e10;
my $outputFormat = "bed"; # bed, detail

my @ARGV0 = @ARGV;

GetOptions ("s"=>\$sameStrand, 
		#"step:i"=>\$step,
		"nms|no-match-score:f"=>\$noMatchScore,
		"of:s"=>\$outputFormat,
		"v"=>\$verbose);


if (@ARGV != 3)
{
	print "extract peak of wiggle values in specified region\n";
	print "Usage: $prog [options] <region.bed> <in.wig> <out.bed>\n";
	print " <in.wig>                 : use - to use stdin; gz file is accepted\n";
	print " -s                       : same strand required\n";
	print " --no-match-score [float] : the score assigned when no match is found ($noMatchScore)\n";
	print " -of             [string] : output format ([bed]|detail)\n"; 
	print " -v                       : verbose\n";
	exit (1);
}



print "CMD=$prog ", join(" ", @ARGV0), "\n" if $verbose;

my ($regionBedFile, $inWiggleFile, $outBedFile) = @ARGV;

print "load bed regions from $regionBedFile ...\n" if $verbose;

my $regions = readBedFile ($regionBedFile, $verbose);

my $n = @$regions;

print "$n regions loaded\n" if $verbose;

my %regionHash;
foreach my $r (@$regions)
{
	my $chrom = $r->{"chrom"};
	Carp::croak "no strand info in ", Dumper ($r), "\n" if $sameStrand && not exists $r->{"strand"};

	my $strand = exists $r->{"strand"} && $sameStrand ? $r->{"strand"} : 'b';

	$r->{"peakScore"} = $noMatchScore;
	$r->{"peakStart"} = $r->{"chromStart"};
	$r->{"peakEnd"} = $r->{"chromEnd"};
	$r->{"halfPeakStart"} = $r->{"chromStart"};
	$r->{"halfPeakEnd"} = $r->{"chromEnd"};

	#for (my $i = 0; $i < $r->{"chromEnd"} - $r->{"chromStart"} + 1; $i++)
	#{
	#	$r->{"score2"}->[$i] = 0;
	#}
	push @{$regionHash{$strand}->{$chrom}}, $r;
}


my $fin;

my $strand = 'b';
my $chrom = "";
my $currRegions = [];
my @currScores;

if ( $inWiggleFile eq "-")
{
    $fin = *STDIN;
}
else
{
    if ($inWiggleFile =~/\.gz$/)
    {
        open ($fin, "gunzip -c $inWiggleFile | ") || Carp::croak "cannot open file $inWiggleFile to read\n";
    }
    else
    {
        open ($fin, "<$inWiggleFile") || Carp::croak "cannot open file $inWiggleFile to read\n";
    }
}

while (my $line = <$fin>)
{
	chomp $line;
	next if $line=~/^\s*$/;
	
	my $firstTrack = 0;
	if ($line=~/^track\s/ || $firstTrack)
	{
		$firstTrack = 0;
		print "a new track found ($line)...\n" if $verbose;

		if (@$currRegions > 0 && @currScores > 0)
		{
			my $nr = @$currRegions;
			my $ns = @currScores;
			print "extracting scores for chrom=$chrom, strand=$strand ($nr regions, $ns scores) ...\n" if $verbose;
			extractPeak ($currRegions, \@currScores); # if @$currRegions > 0 && @currScores > 0;
		}

		$chrom = "";
		$currRegions = [];
		@currScores = ();
		$strand = 'b';

		if ($sameStrand)
	   	{
			if($line=~/\(([\+|\-])\)/)
			{
				$strand = $1;
			}
			else
			{
				Carp::croak "no strand information found in the track\n";
			}
		}
		next if $line =~/^track\s/;
	}
	
	my ($chr, $chromStart, $chromEnd, $score) = split (/\t/, $line);

	Carp::croak "the score in the wiggle can not be less than $noMatchScore\n" if $score <= $noMatchScore;
	#make sure the scores are reasonable

	$chromEnd -= 1;
	
	if ($chr ne $chrom)
	{
		if (@$currRegions > 0 && @currScores > 0)
		{
			my $nr = @$currRegions;
			my $ns = @currScores;
			print "extracting scores for chrom=$chrom, strand=$strand ($nr regions, $ns scores) ...\n" if $verbose;
			extractPeak ($currRegions, \@currScores); # if @$currRegions > 0 && @currScores > 0;
		}

		print "processing $chr ...\n" if $verbose;
		$chrom = $chr;
		$currRegions = [];

		$currRegions = $regionHash{$strand}->{$chrom} if exists $regionHash{$strand} && exists $regionHash{$strand}->{$chrom};
		@currScores = ();
	}
	
	push @currScores, {chrom=>$chr, chromStart=>$chromStart, chromEnd=>$chromEnd, score=>$score};
}

if (@$currRegions > 0 && @currScores > 0)
{
	my $nr = @$currRegions;
	my $ns = @currScores;
	print "extracting scores for chrom=$chrom, strand=$strand ($nr regions, $ns scores) ...\n" if $verbose;
	extractPeak ($currRegions, \@currScores); # if @$currRegions > 0 && @currScores > 0;
}


close ($fin) if $inWiggleFile ne '-';


print "dumping results to $outBedFile ...\n" if $verbose;

if ($outputFormat eq 'bed')
{
	foreach my $r (@$regions)
	{
		#$r->{"chromStart"} = $r->{"peakStart"};
		#$r->{"chromEnd"} = $r->{"peakEnd"};
		$r->{"score"} = $r->{"peakScore"};
	}
	writeBedFile ($regions, $outBedFile);
}
elsif ($outputFormat eq 'detail')
{
	my $fout;
	open ($fout, ">$outBedFile") || Carp::croak "cannot open file $outBedFile to write\n";

	foreach my $r (@$regions)
	{
		print $fout join ("\t", bedToLine ($r), 
				$r->{"peakScore"}, $r->{"peakStart"}, $r->{"peakEnd"} + 1, $r->{"halfPeakStart"}, $r->{"halfPeakEnd"} + 1), "\n";
	}
	close ($fout);
}
else
{
	Carp::croak "wrong output format: $outputFormat\n";
}


sub extractPeak
{
	my ($regions, $scores) = @_;
	my @regionsSorted = sort {$a->{"chromStart"} <=> $b->{"chromStart"}} @$regions;
	
	my $firstScoreIdx = 0;

	print "determine peak height ...\n" if $verbose;

	my $k = 0;
	foreach my $r (@regionsSorted)
	{

		print "$k ...\n" if $verbose && $k % 5000 == 0;
		$k++;

		my $chromStart = $r->{"chromStart"};
		my $chromEnd = $r->{"chromEnd"};
		while ($firstScoreIdx < @$scores && $scores->[$firstScoreIdx]->{"chromEnd"} < $chromStart)
		{
			$firstScoreIdx++;
		}

		#print "firstIdx = $firstScoreIdx\n" if $verbose;

		my $i = $firstScoreIdx;
		while ($i < @$scores && $scores->[$i]->{"chromStart"} <= $chromEnd)
		{
			my $overlapStart = max ($chromStart, $scores->[$i]->{"chromStart"});
			my $overlapEnd = min ($chromEnd, $scores->[$i]->{"chromEnd"});
			
			if ($overlapEnd >= $overlapStart)
			{
				if ($r->{"peakScore"} < $scores->[$i]->{"score"})
				{
					#use the new region
					$r->{"peakScore"} = $scores->[$i]->{"score"};
					$r->{"peakStart"} = $scores->[$i]->{"chromStart"};
					$r->{"peakEnd"} = $scores->[$i]->{"chromEnd"};
				}
				elsif ($r->{"peakScore"} == $scores->[$i]->{"score"})
				{
					#expand the region
					$r->{"peakStart"} = min ($r->{"peakStart"}, $scores->[$i]->{"chromStart"});
					$r->{"peakEnd"} = max ($r->{"peakEnd"}, $scores->[$i]->{"chromEnd"});
				}
			}
			$i++;
		}
	}

	print "determine half peak width ...\n" if $verbose;

	$firstScoreIdx = 0;
	$k = 0;
	foreach my $r (@regionsSorted)
	{

		print "$k ...\n" if $verbose && $k % 5000 == 0;
		$k++;


		my $chromStart = $r->{"chromStart"};
		my $chromEnd = $r->{"chromEnd"};

		next unless $r->{"peakScore"} > $noMatchScore; 


		my $halfPeakStartMatch = 0;

		while ($firstScoreIdx < @$scores && $scores->[$firstScoreIdx]->{"chromEnd"} < $chromStart)
		{
			$firstScoreIdx++;
		}

		#print "firstIdx = $firstScoreIdx\n" if $verbose;

		my $i = $firstScoreIdx;
		while ($i < @$scores && $scores->[$i]->{"chromStart"} <= $chromEnd)
		{
			my $overlapStart = max ($chromStart, $scores->[$i]->{"chromStart"});
			my $overlapEnd = min ($chromEnd, $scores->[$i]->{"chromEnd"});
			
			if ($overlapEnd >= $overlapStart && $scores->[$i]->{"score"} > $r->{"peakScore"} / 2)
			{
				if ($halfPeakStartMatch == 0)
				{
					$r->{"halfPeakStart"} = $scores->[$i]->{"chromStart"};
				    $halfPeakStartMatch = 1;
				}	
				$r->{"halfPeakEnd"} = $scores->[$i]->{"chromEnd"};
			}
			$i++;
		}
	}
}


