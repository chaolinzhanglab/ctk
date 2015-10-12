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
	print "Usage: $prog [options] <in.fa> <out.fa>\n";
	print "options:\n";
	print " --tmp-dir   [string]: temp dir ($tmpDir)\n";
	print " -v                  : verbose\n";
	exit (1);
}

Carp::croak "$tmpDir already exists\n" if -d $tmpDir;

system ("mkdir $tmpDir");


my ($inFastaFile, $outFastaFile) = @ARGV;

my $tmpHeaderFile = "$tmpDir/header";
my $tmpSeqFile = "$tmpDir/seq";

my $cmd = "grep \\> $inFastaFile | cut -d \" \" -f 1 > $tmpHeaderFile";
print $cmd, "\n" if $verbose;

system ($cmd);

$cmd = "grep -v \\> $inFastaFile > $tmpSeqFile";
print $cmd, "\n" if $verbose;
system ($cmd);

$cmd = "paste $tmpHeaderFile $tmpSeqFile | sort -k 2 | uniq -f 1 -c | awk '{print \$2\"#\"\$1\"\\n\"\$3}' > $outFastaFile";
print $cmd, "\n" if $verbose;
system ($cmd);
system ("rm -rf $tmpDir");







