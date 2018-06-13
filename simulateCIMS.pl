#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Carp;

use Bed;

my $prog = basename ($0);

my $verbose = 0;
my $mutationSize = 1;
my $numSites = 0;
my $CIMSPositionFile = "";
my $rand_seed = 0;

GetOptions ("n:i"=>\$numSites,
		'w:i'=>\$mutationSize,
		'p:s'=>\$CIMSPositionFile,
		'srand:f'=>\$rand_seed,
		"v"=>\$verbose);



if (@ARGV != 2)
{
	print "simulate CLIP mismatch sites\n";
	print "Usage: $prog [options] <in.bed> <out.bed>\n";
	print " -n  [int]     : number of mismatches\n";
	print " -w  [int]     : mutation size ($mutationSize)\n";
	print " -p [file]     : specify the position of mismatch sites relative to read start\n";
	print " -srand [float]: seed of random number (0=no fixed seeds)\n";
	print " -v            : verbose\n";
	exit (0);
}


my ($inBedFile, $outBedFile) = @ARGV;


Carp::croak "-n or -p must be specified\n" if $numSites == 0 && $CIMSPositionFile eq '';


my $numTags = `grep -v \"track=\" $inBedFile | wc -l`;
chomp $numTags;
$numTags =~/^\s*(\d+)/;
$numTags = $1;

print "$numTags tags detected in $inBedFile\n" if $verbose;



my @CIMSPositions;
if (-f $CIMSPositionFile)
{
	print "loading CIMS positions ...\n" if $verbose;
	my $fin;
	open ($fin, "<$CIMSPositionFile");
	while (my $line=<$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		Carp::croak "position < 0?\n" if $line < 0;
		push @CIMSPositions, $line;
	}
	$numSites = @CIMSPositions;
}
elsif ($CIMSPositionFile ne '')
{
	Carp::croak "can not read $CIMSPositionFile\n";
}

#record how many time each tag should be sampled and at which position
#if uniform sampling is used, hold the position undetermined

if ($rand_seed != 0)
{
	srand ($rand_seed);
}
my %idxHash;
for (my $i = 0; $i < $numSites; $i++)
{
	my $idx = int(rand($numTags)); #index of tags to sample
	if (@CIMSPositions == $numSites)
	{
		push @{$idxHash{$idx}}, $CIMSPositions[$i];
	}
	else
	{
		push @{$idxHash{$idx}}, -1;
	}
}

print "permute mutations ...\n" if $verbose;
my $fin;
my $fout;

open ($fin, "<$inBedFile") || Carp::croak "cannot open file $inBedFile to read\n";
open ($fout, ">$outBedFile") || Carp::croak "cannot open file $outBedFile to write\n";

my $idx = -1;

my $n = 0;
while (my $line=<$fin>)
{
	chomp $line;
	next if $line=~/^track/;
	next if $line=~/^\s*$/;
	my $t = lineToBed ($line);
	$idx++;

	print "$idx ...\n" if $verbose && $idx % 50000 == 0;


	next unless exists $idxHash{$idx};

	my $positions = $idxHash{$idx};
	foreach my $pos (@$positions)
	{
		my $chromStart = $t->{'chromStart'};
		my $chromEnd = $t->{'chromEnd'};
		if ($pos >= 0)
		{
			$chromStart = $t->{'strand'} eq '+' ? $chromStart + $pos : $chromEnd - $pos;

			if ($chromStart < 0)
			{
				#cz fix 06/12/2018
				#in some very rare cases, the mutation go outside the range
				#as a quick and dirty fix, we do uniform sampling for these
				
				#uniform sampling
            	$chromStart = $t->{'chromStart'} + 1; #do not allow changes in the first and last position
            	$chromEnd = $t->{'chromEnd'} - 1;
            	my $tagLen = $chromEnd - $chromStart + 1;
            	$pos = int(rand($tagLen - $mutationSize + 1));
            	$chromStart += $pos;
			}
		}
		else
		{
			#uniform sampling
			$chromStart = $t->{'chromStart'} + 1; #do not allow changes in the first and last position
			$chromEnd = $t->{'chromEnd'} - 1;
			my $tagLen = $chromEnd - $chromStart + 1;
			$pos = int(rand($tagLen - $mutationSize + 1));
			$chromStart += $pos;
		}
		print $fout join ("\t", $t->{'chrom'}, $chromStart, $chromStart + $mutationSize, "CIM_". $n, 0, $t->{'strand'}), "\n";
		$n++;
	}
}

Carp::croak "the total number of mutation sites is not $numSites as expected\n" unless $n == $numSites;

close ($fout);
close ($fin);

