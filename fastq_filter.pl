#!/usr/bin/perl -w
#

use strict;
use File::Basename;
use Getopt::Long;
use Carp;
use Data::Dumper;

use Common;


my $prog = basename($0);
my $verbose = 0;
my $fastqFormat = "sanger"; #64 offset, or 'sanger' for 33 offset


my $filterStr = "";
my @filters;
my $maxN = -1;
my $outputFormat = "fasta";

my $index = "";
my $indexStart = 0; #0-based
my $indexStr = "";

#my $uniqFlag = 0;

my @ARGV0 = @ARGV;

GetOptions (
			'if:s'=>\$fastqFormat,
			'f:s'=>\$filterStr,
			'index:s'=>\$index,
			#'uniq'=>\$uniqFlag,
			'maxN:i'=>\$maxN,
			'of:s'=>\$outputFormat,
			'v|verbose'=>\$verbose
);

if (@ARGV != 2)
{
	print "Convert solexa raw data to fasta file\n";
	print "Usage: $prog [options] <in.fq> <out.fa or output.fq>\n";
	print "You can also use - to specify STDIN for input or STDOUT for output\n\n";
	print "OPTIONS:\n";
	print " -v              : verbose\n";
	print " -if    [string] : fastq format (solexa|[sanger])\n";
	print " -index [string] : index position and sequence (e.g. 0:CAGT)\n";
	print " -f     [string] : quality score filter string \n";
	print "                   (method1:start1-end1:score1,method2:start2-end2:score2, zero-based)\n";
	print " -maxN  [int]    : max number of N in sequence (default off)\n";
	print " -of    [string] : output format ([fasta]|fastq)\n";
	print "\n";
	exit (1);
}


my ($in, $out) = @ARGV;

my $offset;
if ($fastqFormat eq 'solexa')
{
	$offset = 64;
}
elsif ($fastqFormat eq 'sanger')
{
	$offset = 33;
}
else
{
	Carp::croak "incorrect fastq format:$fastqFormat\n";
}


print STDERR "CMD = $prog ", join (" ", @ARGV0), "\n" if $verbose;

#filter by quality scores
my @cols = split (",", $filterStr);
foreach my $f (@cols)
{
	$f=~/^(\S*?):(\d+)\-(\d+):(\S+)/;
	
	print STDERR "f=$f\n" if $verbose;
	
	my $ff = {start=>$2, end=>$3, score=>$4};
	my $subName = $1;
	
	if ($subName eq 'mean')
	{
		$ff->{"sub"} = \&mean;
	}
	elsif ($subName eq 'min')
	{
		$ff->{"sub"} = \&min;
	}
	else
	{
		Carp::croak "incorrect subroutine in filter: $subName\n";
	}
	push @filters, $ff;
}

my $nfilters = @filters;

print STDERR "$nfilters filters detected\n" if $verbose;

#filter by index
if ($index ne '')
{
	($indexStart, $indexStr) = split (":", $index);
	$indexStr = uc ($indexStr);

	if ($indexStr eq '' || $indexStr=~/[^ACGT]/)
	{
		Carp::croak "incorrect format of index:$index\n";
	}
}


#my ($seqId, $seq, $score) = ("", "", "");

my ($fin, $fout);

if ( $in eq "-")
{
    $fin = *STDIN;
}
else
{
	if ($in =~/\.gz$/)
	{
		open ($fin, "gunzip -c $in | ") || Carp::croak "cannot open file $in to read\n";
	}
	else
	{
    	open ($fin, "<$in") || Carp::croak "cannot open file $in to read\n";
	}
}
if ( $out eq "-")
{
     $fout = *STDOUT;
}
else
{
    open ($fout, ">$out") || Carp::croak "cannot open file $out to write\n";
}

my $iter = 0;

my $indexLen = $indexStr ne '' ? length($indexStr) : 0;

while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;
	
	#expect only the id line
	Carp::croak "the id line not found\n" unless $line =~/^\@(.*?)$/;
	my $seqId = $1;
	$seqId =~s/\s.*$//g;
	
	#print $seqId, "\n";
	
	#read the next three lines, we do not expect blank lines in between
	my $seq = <$fin>;
	chomp $seq;

	my $tmp = <$fin>;
	my $qualityScoreLine = <$fin>;
	chomp $qualityScoreLine;

	#print $qualityScoreLine, "\n";
	Carp::croak "data inconsistency found\n" unless $seq && $qualityScoreLine;

	if ($index ne '')
	{
		next unless uc (substr ($seq, $indexStart, $indexLen)) eq $indexStr;
	}


	#parse the quality scores
	my @score;
	if ($qualityScoreLine =~/^\s*\-?\d+\s/)
	{
		#numerical format
		$qualityScoreLine =~s/^\s+//g;#remove the leading space if any
		@score = split (/\s+/, $qualityScoreLine);
		#print join ("\t", @score), "\n";
	}
	else
	{
		#ascii format
		@score = split (//, $qualityScoreLine);
		@score = map {ord ($_) - $offset} @score;
		#print join ("\t", @score), "\n";
	}

	my $pass = 1;
	foreach my $f (@filters)
	{
		#print Dumper ($f), "\n";
		my @s = @score[$f->{"start"} .. $f->{"end"}];
		my $ret = $f->{"sub"}->(\@s);

		#print "ret = $ret, score=", join ("\t", @s), "\n";
		if ($ret < $f->{"score"})
		{
			$pass = 0;
			last;
		}
	}

	if ($maxN >= 0)
	{
		my $goodBase = ($seq=~tr/ACGTacgt//);
		my $n = length($seq) - $goodBase;
		$pass = 0 if $n > $maxN;
	}

	if ($pass)
	{
		#$seqId =~s/\/\d+$//g;
		if ($outputFormat eq 'fasta')
		{
			print $fout ">$seqId\n";
			print $fout $seq, "\n";
		}
		elsif ($outputFormat eq 'fastq')
		{
			
			print $fout "\@$seqId\n";
			print $fout $seq, "\n";
			print $fout "+\n";
			print $fout $qualityScoreLine, "\n";
		}
		else
		{
			Carp::croak "unknown output format: $outputFormat\n";
		}
	}

	#($seqId, $seq, $score) = ("", "", "");
	$iter++;
		
	print STDERR "$iter ...\n" if $verbose && $iter % 100000 == 0;
}

close ($fin) if $in ne '-';
close ($fout);

