#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Carp;
use File::Basename;

use MyConfig;

my $linkerLen = 5;
#my $seqLen = 36;
my $verbose = 0;
my $prog = basename ($0);
my $format = "fasta"; # or 'fastq'
my $tmpDir = getDefaultCache ($prog);
my $keepCache = 0;

my $linkerStartWith = "";
my $linkerEndWith = "";

my @ARGV0 = @ARGV;

GetOptions ("len:i"=>\$linkerLen,
		#"seq:i"=>\$seqLen,
		"format:s"=>\$format,
		"barcode-start-with:s"=>\$linkerStartWith,
		"barcode-end-with:s"=>\$linkerEndWith,
		"tmp-dir:s"=>\$tmpDir,
		'keep-cache'=>\$keepCache,
		"v|verbose"=>\$verbose);

if (@ARGV != 2)
{
	print "strip random linker sequences and attach that in sequence ids\n";
	print "Usage: $prog [options] <file.in> <file.out>\n";
	print " <file.in> : fasta or fastq file, .gz file accepted\n";
	print "options:\n";
	print " -len       [int]             : length of barcode sequences ($linkerLen)\n";
	print " -format    [string]          : input format ([fasta]|fastq)\n";
	print " --barcode-start-with [string]: filter sequences based on the starting nucleotides in the barcode\n";
	print " --barcode-end-with   [string]: filter sequences based on the ending nucleotides in the barcode\n";
	#print " -seq           [int]: sequence length ($seqLen)\n";
	print " --tmp-dir   [string]         : temp dir ($tmpDir)\n";
	print " --keep-cache                 : keep cache when the job is done\n";
	print " -v                           : verbose\n";
	exit (1);
}



print "CMD = $prog ", join (" ", @ARGV0), "\n" if $verbose;


#Carp::croak "$tmpDir already exists\n" if -d $tmpDir;

system ("mkdir $tmpDir");


my $linkerStartFilterLen = length($linkerStartWith);
my $linkerEndFilterLen = length ($linkerEndWith);

Carp::croak "errors in filter patterns\n" if $linkerStartFilterLen > $linkerLen || $linkerEndFilterLen > $linkerLen;

my ($inFile, $outFile) = @ARGV;

if ($format eq 'fasta')
{
	my $tmpHeaderFile = "$tmpDir/header";
	my $tmpSeqFile = "$tmpDir/seq";

	my $cmd = "grep \\> $inFile > $tmpHeaderFile";
	$cmd = "gunzip -c $inFile | grep \\> > $tmpHeaderFile" if $inFile =~/\.gz$/;
		
	print $cmd, "\n" if $verbose;

	system ($cmd);

	$cmd = "grep -v \\> $inFile > $tmpSeqFile";
	$cmd = "gunzip -c $inFile | grep -v \\> > $tmpSeqFile" if $inFile =~/\.gz$/;

	print $cmd, "\n" if $verbose;
	system ($cmd);

	my $tmpSeqLinkerFile = "$tmpDir/linker";
	$cmd = "cut -c 1-$linkerLen $tmpSeqFile > $tmpSeqLinkerFile";
	print $cmd, "\n" if $verbose;
	system ($cmd);

	my $tmpSeqReadFile = "$tmpDir/reads";
	$cmd = "cut -c " . ($linkerLen + 1) . "- $tmpSeqFile > $tmpSeqReadFile";
	print $cmd, "\n" if $verbose;
	system ($cmd);

	#$cmd = "paste $tmpHeaderFile $tmpSeqLinkerFile $tmpSeqReadFile | awk '{print \$1\"#\"\$2\"\\n\"\$3}' > $outFile";#no filter
	$cmd = "paste $tmpHeaderFile $tmpSeqLinkerFile $tmpSeqReadFile | awk '{if(1";
	$cmd .= " && substr(\$2, 1, $linkerStartFilterLen) == \"$linkerStartWith\"" if $linkerStartWith ne '';
	$cmd .= " && substr(\$2, ".($linkerLen - $linkerEndFilterLen + 1) . ", $linkerEndFilterLen) == \"$linkerEndWith\"" if $linkerEndWith ne '';
	$cmd .= ") {print \$1\"#\"\$2\"\\n\"\$3}}' > $outFile";

	#Carp::croak $cmd, "\n";
	print $cmd, "\n" if $verbose;
	system ($cmd);
}
elsif ($format eq 'fastq')
{
	my $tmpIdFile = "$tmpDir/id";
	
	#my $cmd = "grep  -v \"^[^\@\+]\" $inFile > $tmpIdFile";
	#$cmd = "gunzip -c $inFile | grep  -v \"^[^\@\+]\" > $tmpIdFile" if $inFile =~/\.gz$/;
	
	#extract two id lines
	my $cmd = "awk '{if(NR%4==1 || NR%4==3) {print \$1}}' $inFile > $tmpIdFile";
	$cmd = "gunzip -c $inFile | awk '{if(NR%4==1 || NR%4==3) {print \$1}}' > $tmpIdFile" if $inFile =~/\.gz$/;

	print $cmd, "\n" if $verbose;
	system ($cmd);

	#extract sequences
	my $tmpSeqFile = "$tmpDir/seq";
	#$cmd = "grep  \"^[^\@\+]\" $inFile | awk 'BEGIN {i=0} {if (i%2==0) {print \$0};i=i+1}' > $tmpSeqFile";
	#$cmd = "gunzip -c $inFile | grep  \"^[^\@\+]\" | awk 'BEGIN {i=0} {if (i%2==0) {print \$0};i=i+1}' > $tmpSeqFile" if $inFile =~/\.gz$/;
	$cmd = "awk '{if(NR%4==2) {print \$1}}' $inFile > $tmpSeqFile";
	$cmd = "gunzip -c $inFile | awk '{if(NR%4==2) {print \$1}}' > $tmpSeqFile" if $inFile =~/\.gz$/;

	print $cmd, "\n" if $verbose;
	system ($cmd);

	#extract linker sequences, note that it is duplicated for the two id lines
	my $tmpSeqLinkerFile = "$tmpDir/linker";
	$cmd = "cut -c 1-$linkerLen $tmpSeqFile | awk '{print \$0\"\\n\"\$0}'> $tmpSeqLinkerFile";
	print $cmd, "\n" if $verbose;
	system ($cmd);

	my $tmpIdFile2 = "$tmpDir/id2";
	$cmd = "paste -d \"#\" $tmpIdFile $tmpSeqLinkerFile > $tmpIdFile2";
	print $cmd, "\n" if $verbose;
	system ($cmd);


	my $tmpSeqReadFile = "$tmpDir/reads";
	$cmd = "cut -c " . ($linkerLen + 1) . "- $tmpSeqFile > $tmpSeqReadFile";
	print $cmd, "\n" if $verbose;
	system ($cmd);



	my $tmpQualityFile = "$tmpDir/quality";
	#$cmd = "grep  \"^[^\@\+]\" $inFile | awk 'BEGIN {i=0} {if (i%2!=0) {print \$0};i=i+1}' > $tmpQualityFile";
	#$cmd = "gunzip -c $inFile | grep  \"^[^\@\+]\" | awk 'BEGIN {i=0} {if (i%2!=0) {print \$0};i=i+1}' > $tmpQualityFile" if $inFile =~/\.gz$/;
	$cmd = "awk '{if(NR%4==0) {print \$1}}' $inFile > $tmpQualityFile";	
	$cmd = "gunzip -c $inFile | awk '{if(NR%4==0) {print \$1}}' > $tmpQualityFile" if $inFile =~/\.gz$/;

	print $cmd, "\n" if $verbose;
	system ($cmd);

	my $line = `head -n 1 $tmpQualityFile`;
	chomp $line;
	my $tmpReadQualityFile = "$tmpDir/read_quality";

	if ($line=~/\s/)
	{
		#numerical score
		$cmd = "cut -f " . ($linkerLen + 1) . "- $tmpQualityFile > $tmpReadQualityFile";	
	}
	else
	{
		#ascii score
		$cmd = "cut -c " . ($linkerLen + 1) . "- $tmpQualityFile > $tmpReadQualityFile";
	}
	print $cmd, "\n" if $verbose;
	system ($cmd);

	$cmd = "paste -d \"\\n\" $tmpSeqReadFile $tmpReadQualityFile > $tmpSeqFile";
	print $cmd, "\n" if $verbose;
	system ($cmd);

	#print only 1 id line
	$cmd = "paste $tmpIdFile2 $tmpSeqFile | awk '{if(NR%2==1){print \$1\"\\n\"\$2}else{print \"+\\n\"\$2}}'> $outFile";
	
	if ($linkerStartWith ne '' || $linkerEndWith ne '')
	{
		my $linkerLenMinus1 = $linkerLen - 1;
		$cmd = "paste $tmpIdFile2 $tmpSeqFile | awk '{linker=substr(\$1,length(\$1)-$linkerLenMinus1, $linkerLen); if(1";
		$cmd .= " && substr(linker, 1, $linkerStartFilterLen) == \"$linkerStartWith\"" if $linkerStartWith ne '';
		$cmd .= " && substr(linker, ".($linkerLen - $linkerEndFilterLen + 1) . ", $linkerEndFilterLen) == \"$linkerEndWith\"" if $linkerEndWith ne '';
		#$cmd .= ") {print \$1\"\\n\"\$2}}' > $outFile";
		$cmd .= ") {if(NR%2==1){print \$1\"\\n\"\$2}else{print \"+\\n\"\$2}}}' > $outFile";
	}
	#Carp::croak $cmd, "\n";
	print $cmd, "\n" if $verbose;
	system ($cmd);
}
else
{
	Carp::croak "unknown format: $format\n";
}
system ("rm -rf $tmpDir") unless $keepCache;






