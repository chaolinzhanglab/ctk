#!/usr/bin/perl -w
#

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Carp;


use Bed;

my $prog = basename($0);

my $maskRepeats = 0;
my $neighbor = ''; # or 'up' or 'down' or 'r=100'
my ($left, $right) = (0, 0);
my $chrLen = {};
my $chrLenFile = "";

my $verbose = 0;


GetOptions (
			'n|neighbor:s'=>\$neighbor,
			'l|left:i'=>\$left,
			'r|right:i'=>\$right,
			'chrLen:s'=>\$chrLenFile,
			'v|verbose'=>\$verbose
);

if (@ARGV != 2)
{
	print "Extract sequences specied by a bed file\n";
	print "Usage: $prog [options] <in.bed> <out.bed>\n";
	print " use \"-\" for stdin or stdout\n";
	print "OPTIONS:\n";
	print " -n : get neighbor region relative to <up|down|r=100>\n";
	print " -l : extension on the left, with sign <$left>\n";
	print " -r : extension on the right, with sign <$right>\n";
	print " -chrLen     [string]: chromosome length file\n";
	print " -v                  : verbose\n";
	print "\n";
	exit (1);
}


my ($in, $out) = @ARGV;

if ($chrLenFile)
{
	my $fin;
	open ($fin, "<$chrLenFile") || Carp::croak "can not open file $chrLenFile to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		my ($chrom, $len) = split ("\t", $line);
		$chrLen->{$chrom} = $len;
	}
	close ($fin);
}


#print "loading regions from $in ...\n";
#my $regions = readBedFile($in);

my $fin;

if ($in eq '-')
{
	$fin = *STDIN;
}
else
{
	if ($in =~/\.gz$/)
	{
		open ($fin, "gunzip -c $in |") || Carp::croak "cannot open file $in to read\n";
	}
	else
	{
		open ($fin, "<$in") || Carp::croak "cannot open file $in to read\n";
	}
}

my $foutBed;
{
	Carp::croak "$out already exists\n" if -f $out && $out ne "-";

	if ($out eq '-')
	{
		$foutBed = *STDOUT;
	}
	else
	{
		open ($foutBed, ">$out") || Carp::croak "can not open file $out to write\n";
	}
}

my $i = 0;
while (my $line = <$fin>)
{
	chomp $line;

	next if $line =~/^\s*$/;

	print STDERR "$i ...\n" if $verbose && $i % 100000 == 0;
	$i++;

	my $r = lineToBed ($line);

	my $chrom = $r->{"chrom"};
	my $chromStart = $r->{"chromStart"};
	my $chromEnd = $r->{"chromEnd"};
	my $name = $r->{"name"} || "$chrom:$chromStart-" . ($chromEnd+1);
	my $strand = '+';
	if (exists $r->{"strand"})
	{
		$strand = '-' if $r->{"strand"} eq '-';
	}

	#modified to allow extension
	if ($neighbor eq '')
	{
		if ($strand eq '+')
		{
			$chromStart += $left;
			$chromEnd += $right;
		}
		else
		{
			$chromStart -= $right;
			$chromEnd -= $left;
		}
		if ($chromEnd < $chromStart)
		{
			print STDERR "$chrom:$chromStart-$chromEnd is too small for truncation\n" if $verbose;
			next;
		}
	}
	elsif ($neighbor eq 'down' && $right - $left >= 0)
	{
		if ($strand eq '+')
		{
			$chromStart = $chromEnd + $left;
			$chromEnd = $chromEnd + $right;
		}
		else
		{
			$chromEnd = $chromStart - $left;
			$chromStart = $chromStart - $right;
		}
	}
	elsif ($neighbor eq 'up' && $right - $left >= 0)
	{
		if ($strand eq '+')
		{
			$chromEnd = $chromStart + $right;
			$chromStart = $chromStart + $left;
		}
		else
		{
			$chromStart = $chromEnd - $right;
			$chromEnd = $chromEnd - $left;
		}
	}
	elsif ($neighbor =~/r=(\d+)$/)
	{
		my $sampleLen = $1;
		my $seqLen = $chromEnd - $chromStart + 1;
		warn "$name (length=$seqLen < $sampleLen) is too small\n" if $seqLen < $sampleLen;
		
		my $offset = int(rand($seqLen - $sampleLen));
		$chromStart += $offset;
		$chromEnd = $chromStart + $sampleLen - 1;
	}
	else
	{
		Carp::croak "which region do you want?\n";
	}

	if (exists $chrLen->{$chrom})
	{
		if ($chromStart < 0)
		{
			$chromStart = 0;
			$chromEnd = $chromStart if $chromEnd < $chromStart;
		}
		
		if ($chromEnd >= $chrLen->{$chrom})
		{
			$chromEnd = $chrLen->{$chrom} -1;
			$chromStart = $chromEnd if $chromStart > $chromEnd;
		}
	}

	my $r2 = {chrom=>$r->{"chrom"},
			chromStart=> $chromStart,
			chromEnd=> $chromEnd,
			name=>$name,
			strand=> ($strand eq '+')? '+' : '-'
	};
	$r2->{"score"} = (exists $r->{"score"})? $r->{"score"} : 0;

	print $foutBed bedToLine ($r2), "\n" 
	if $chromEnd >= $chromStart; #fix March 14, 2013
}

close ($foutBed) if $out ne '-';
close ($fin) if $in ne '-';


