#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Carp;
use File::Basename;

use MyConfig;

my $verbose = 0;
my $prog = basename ($0);
my $tmpDir = getDefaultCache ($prog);


GetOptions (
		"tmp-dir:s"=>\$tmpDir,
		"v|verbose"=>\$verbose);

if (@ARGV != 2)
{
	print "collapse exact duplicate sequences\n";
	print "Usage: $prog [options] <in.fq> <out.fq>\n";
	print " use - for stdout\n";
	print "options:\n";
	print " --tmp-dir   [string]: temp dir ($tmpDir)\n";
	print " -v                  : verbose\n";
	exit (1);
}

Carp::croak "$tmpDir already exists\n" if -d $tmpDir;

system ("mkdir $tmpDir");


my ($inFastqFile, $outFastqFile) = @ARGV;

my $tmpIDFile = "$tmpDir/id";
my $tmpSeqFile = "$tmpDir/seq";
my $tmpQualFile = "$tmpDir/qual";

if ($inFastqFile =~/\.gz$/)
{
	my $cmd = "gunzip -c $inFastqFile | awk '{if(NR%4==1) {print \$1}}' > $tmpIDFile";
	print $cmd, "\n" if $verbose;
	system ($cmd);

	$cmd = "gunzip -c $inFastqFile | awk '{if(NR%4==2) {print \$0}}' > $tmpSeqFile";
	print $cmd, "\n" if $verbose;
	system ($cmd);

	$cmd = "gunzip -c $inFastqFile | awk '{if(NR%4==0) {print \$0}}' > $tmpQualFile";
	print $cmd, "\n" if $verbose;
	system ($cmd);
}
else
{
	my $cmd = "awk '{if(NR%4==1) {print \$1}}' $inFastqFile > $tmpIDFile";
	print $cmd, "\n" if $verbose;
	system ($cmd);

	$cmd = "awk '{if(NR%4==2) {print \$0}}' $inFastqFile > $tmpSeqFile";
	print $cmd, "\n" if $verbose;
	system ($cmd);

	$cmd = "awk '{if(NR%4==0) {print \$0}}' $inFastqFile > $tmpQualFile";
	print $cmd, "\n" if $verbose;
	system ($cmd);
}

my $cmd = "paste $tmpIDFile $tmpQualFile $tmpSeqFile| sort -T $tmpDir -k 3 | uniq -f 2 -c | awk '{print \$2\"#\"\$1\"\\n\"\$4\"\\n+\\n\"\$3}' > $outFastqFile";
$cmd = "paste $tmpIDFile $tmpQualFile $tmpSeqFile| sort -T $tmpDir -k 3 | uniq -f 2 -c | awk '{print \$2\"#\"\$1\"\\n\"\$4\"\\n+\\n\"\$3}'" if $outFastqFile eq '-';

print $cmd, "\n" if $verbose;
system ($cmd);
system ("rm -rf $tmpDir");







