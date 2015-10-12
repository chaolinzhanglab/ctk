#!/usr/bin/perl -w


use strict;
use warnings;

use File::Basename;
use Data::Dumper;
use Getopt::Long;
use Carp;

my $prog = basename($0);

my $mapFileList = "";
my $mismatchFile = "";
my $randomizeMismatch = 0;
my $inDelInScore = 0;

#my $dir = "";
#my $mapFile = ""; 
#my $readLen = 26; #length of mapped  reads
#my $maxLen = 0; #maximize length <=2 mismatches instead of minimize mismatches
my $verbose = 0;

GetOptions (#'list:s'=>\$mapFileList,
			#'d:s'=>\$dir,
			#'map:s'=>\$mapFile,
			#'l|read-length:i'=>\$readLen,
			#'max-len'=>\$maxLen,
			'mismatch-file:s'=>\$mismatchFile,
			'randomize-mismatch'=>\$randomizeMismatch,
			'indel-in-score'=>\$inDelInScore,
			'v|verbose'=>\$verbose
			);

if (@ARGV != 2)
{
	print "convert and combine novo align output to bed format\n";
	print "Usage: $prog [options] <in.novoalign> <out.bed>\n";
	print " --mismatch-file [string]: mismatch-file\n";
	print " --indel-in-score        : report the number of mismatches including indels in the score column\n";
	#print " --randomize-mismatch    :\n";
	#print " -d    [string]: directory where map files are located\n";
	#print " -list [file]  : list of mapping file names (file/read size pairs)\n";
	#print " -max-len      : max length <=2 mismatches\n";
	#print " -map  [file]  : map file\n";
	#print " -l    [int]   : read length ($readLen)\n";
	print " -v                      : verbose (on|[off])\n";
	exit (1);
}


my ($mapFile, $outBedFile) = @ARGV;


#get the number of sequences in the file
my $seqNo = `wc -l $mapFile`;
$seqNo =~/^\s*(\d+)\s/;
$seqNo = $1;

print "$seqNo reads to be combined\n" if $verbose;

my $fin;
open ($fin, "<$mapFile") || Carp::croak "cannot open file $mapFile to read\n";

my $fout;
open ($fout, ">$outBedFile") || Carp::croak "can not open file $outBedFile to write\n";

my $foutMismatch;
if ($mismatchFile ne '')
{
	open ($foutMismatch, ">$mismatchFile") || Carp::croak "can not open file $mismatchFile to write\n";
}

my $i = 0;
while (my $line = <$fin>)
{
	print "$i ...\n" if $verbose && $i - int($i / 100000) * 100000 == 0;
	$i++;

	chomp $line;
	next if $line=~/^\t*$/;
	next if $line=~/^#/;

	my @cols = split (/\t/, $line);

	#my $n = @cols;
	#print "n=$n\n";
	#Carp::croak join ("\n", @cols), "\n";

	my $name = substr($cols[0], 1);
	
	my $seq = $cols[2];
	my $mapFlag = $cols[4];
	
	next unless $mapFlag eq 'U';

	my $chrom = substr($cols[7], 1);
	my $chromStart = $cols[8] -1;
	my $strand = $cols[9] eq 'F' ? '+' : '-';
	my $mapDetails = @cols < 14 ? '' : $cols[13];

	my $mapLen0 = length ($seq);
	my $mapLen = $mapLen0;
	
	my @items = split (/\s/, $mapDetails);
	my $mismatch = 0;
	if (@items > 0)
	{
		foreach my $item (@items)
		{
			my $matchStart = 1; #always 1 for novoalign
			if ($item=~/^(\d+)(\S)\>(\S)$/)
			{
				#substitution
				my $itemStartOnRead = $randomizeMismatch ? int(rand($mapLen)) : ($1 - 1);
				
				my $itemStart = $chromStart + $itemStartOnRead;
				my $refBase = $2;
				my $mismatchBase = $3;
				if ($itemStartOnRead > 0 && $itemStartOnRead < $mapLen - 1)
				{
					$mismatch++;
					print $foutMismatch join ("\t", $chrom, $itemStart, $itemStart+1, $name, $itemStartOnRead, $strand, $itemStartOnRead, $refBase, ">", $mismatchBase, $matchStart), "\n" if -f $mismatchFile;
				}
			}
			elsif ($item=~/^(\d+)\-(\S)$/) #deletion in reads
			{
				$mapLen++;
				my $refBase = $2;
				my $itemStartOnRead = $randomizeMismatch ? int(rand($mapLen)) : ($1 - 1);
				
				my $itemStart = $chromStart + $itemStartOnRead;
				if ($itemStartOnRead > 0 && $itemStartOnRead < $mapLen - 1)
				{
					$mismatch++ if $inDelInScore;
					print $foutMismatch join ("\t", $chrom, $itemStart, $itemStart+1, $name, $itemStartOnRead, $strand, $itemStartOnRead, $refBase, "-", ".", $matchStart), "\n" if -f $mismatchFile;
				}
			}
			elsif ($item=~/^(\d+)\+(\S+)$/)
			{
				#insertions
				my $itemStartOnRead = $randomizeMismatch ? int(rand($mapLen)) : ($1 - 1);#insertion before this position
				my $refBase = ".";

				my $itemStart = $chromStart + $itemStartOnRead;
				my $mismatchBase = $2;
				$mapLen -= length ($mismatchBase);
				if ($itemStartOnRead > 0 && $itemStartOnRead < $mapLen - 1)
				{
					$mismatch++ if $inDelInScore;
					print $foutMismatch join ("\t", $chrom, $itemStart, $itemStart+1, $name, $itemStartOnRead, $strand, $itemStartOnRead, $refBase, "+", $mismatchBase, $matchStart), "\n" if -f $mismatchFile;
				}
			}
			else
			{
				Carp::croak "pattern not handled properly: $item in $name\n";
			}
		}
	}

	my $chromEnd = $chromStart + $mapLen;
	#my $mismatch = @items;
	#$mismatch -= $mismatchAtEnd;
	print $fout join ("\t", $chrom, $chromStart, $chromEnd, $name, $mismatch, $strand), "\n" if $chromStart >= 0;
}

close ($fin);
close ($fout);
close ($foutMismatch) if $mismatchFile ne '';


