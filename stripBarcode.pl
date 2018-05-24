#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Carp;
use File::Basename;

use MyConfig;

my $barcodeLen = 5;
#my $seqLen = 36;
my $verbose = 0;
my $prog = basename ($0);
my $format = "fasta"; # or 'fastq'

my $barcodeStartWith = "";
my $barcodeEndWith = "";

my @ARGV0 = @ARGV;

GetOptions ("len:i"=>\$barcodeLen,
		#"seq:i"=>\$seqLen,
		"format:s"=>\$format,
		"barcode-start-with:s"=>\$barcodeStartWith,
		"barcode-end-with:s"=>\$barcodeEndWith,
		"v|verbose"=>\$verbose);

if (@ARGV != 2)
{
	print "strip random linker sequences and attach that in sequence ids\n";
	print "Usage: $prog [options] <file.in> <file.out>\n";
	print " <file.in> : fasta or fastq file, .gz file accepted\n";
	print " <file.out>: use - for stdout\n";
	print " Note: numeric quality scores are not supported anymore\n";
	print "options:\n";
	print " -len       [int]             : length of barcode sequences ($barcodeLen)\n";
	print " -format    [string]          : input format ([fasta]|fastq)\n";
	print " --barcode-start-with [string]: filter sequences based on the starting nucleotides in the barcode\n";
	print " --barcode-end-with   [string]: filter sequences based on the ending nucleotides in the barcode\n";
	print " -v                           : verbose\n";
	exit (1);
}


my ($inFile, $outFile) = @ARGV;

my $msgio = $outFile eq '-' ? *STDERR :  *STDOUT;

print $msgio "CMD = $prog ", join (" ", @ARGV0), "\n" if $verbose;


#Carp::croak "$tmpDir already exists\n" if -d $tmpDir;

$barcodeStartWith = uc ($barcodeStartWith) if $barcodeStartWith ne '';
$barcodeEndWith = uc ($barcodeEndWith) if $barcodeEndWith ne '';

my $barcodeStartFilterLen = length($barcodeStartWith);
my $barcodeEndFilterLen = length ($barcodeEndWith);

Carp::croak "errors in filter patterns\n" if $barcodeStartFilterLen > $barcodeLen || $barcodeEndFilterLen > $barcodeLen;


my ($fin, $fout);

if ($inFile eq '-')
{
	$fin = *STDIN;
}
elsif ($inFile =~/\.gz$/)
{
 	open ($fin, "gunzip -c $inFile | ")||Carp::croak "cannot open file $inFile to read\n";
}
else
{
	open ($fin, "<$inFile") || Carp::croak "cannot open file $inFile\n";	
}

if ($outFile eq '-')
{
	$fout = *STDOUT;
}
else
{
	open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";
}

my $iter = 0;

while (my $line = <$fin>)
{
	chomp $line;
	next if $line=~/^\s*$/;
	my ($id, $seq, $qual);
	
	$id = substr ($line, 1);
	if ($format eq 'fasta')
	{
		$seq = <$fin>; chomp $seq;
	}
	else
	{
		#fastq
		$seq = <$fin>; chomp $seq;
		<$fin>;
		$qual = <$fin>; chomp $qual;
		Carp::croak "cannot have space in quality score line\n" if $qual=~/^\s/;
		Carp::croak "seq and qual do not have the same size\n" if length ($seq) != length ($qual);
	}
	
	print $msgio "$iter ...\n" if $verbose && $iter % 100000 == 0;

	$iter++;

	next unless length($seq) > $barcodeLen;	

	my $barcode = substr ($seq, 0, $barcodeLen);
	$id = "$id#$barcode";

	$seq = substr ($seq, $barcodeLen);

	$qual = substr ($qual, $barcodeLen) if $format eq 'fastq';

	if ($barcodeStartWith ne '')
	{
		next unless uc (substr($barcode, 0, $barcodeStartFilterLen)) eq $barcodeStartWith;
	}
	
	if ($barcodeEndWith ne '')
	{
		next unless uc (substr($barcode, $barcodeLen - $barcodeEndFilterLen, $barcodeEndFilterLen)) eq $barcodeEndWith;
	}
	
	if ($format eq 'fasta')
	{
		print $fout ">$id\n";
		print $fout $seq, "\n";
	}
	else
	{
		print $fout "@", $id, "\n";
		print $fout $seq, "\n";
		print $fout "+\n";
		print $fout $qual, "\n";
	}
}

close ($fin) unless $inFile eq '-';
close ($fout) unless $outFile eq '-';






