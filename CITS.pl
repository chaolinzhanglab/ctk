#!/usr/bin/perl -w
#
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Carp;
use Data::Dumper;

use Common;
use MyConfig;

my $prog = basename ($0);
my $cmdDir = dirname ($0);


my $bigFile = "";        #if yes, we need to use cache
#my $minBlockSize = 2000000;
my $pvalueThreshold = 0.01;
my $multiTestCorrection = "";
my $maxGap = -1;
my $prefix = "CITS";

my $cache = getDefaultCache ($prog);
my $keepCache = "";
my $verbose = "";


GetOptions (
		'big'=>\$bigFile,
        'p:f'=>\$pvalueThreshold,
        'multi-test'=>\$multiTestCorrection,
        'gap:i'=>\$maxGap,
        'c|cache:s'=>\$cache,
        'keep-cache'=>\$keepCache,
        'v|verbose'=>\$verbose
);

if (@ARGV != 3)
{
	print "identify significant CITS\n";
	print "usage: $prog [options] <uniq.tag.bed> <uniq.mutation.type.bed> <out.CITS.bed>\n";
	print " -big          : the input tag file is big\n";
	print " -p   [double] : p-value threshold ($pvalueThreshold)\n";
	print " --multi-test  : perform bonferroni multiple test correction\n";
	print " --gap   [int] : max gap used to cluster CITS ($maxGap, -1=no cluster)\n";
	print " -c   [string] : cache dir ($cache)\n";
	print " -v            : verbose\n";
	exit (1);
}

my ($uniqTagBedFile, $uniqMutationBedFile, $outCITSBedFile)=@ARGV;

my $ret = system ("mkdir $cache");
Carp::croak "cannot mkdir $cache\n" if $ret != 0;

##
print "remove tags with potential CITS ...\n" if $verbose;

my $uniqTagCleanBedFile = "$cache/tag.clean.bed";

my $cmd = "perl $cmdDir/removeRow.pl -q 3 -f 3 $uniqTagBedFile $uniqMutationBedFile > $uniqTagCleanBedFile";
print $cmd, "\n" if $verbose;

$ret = system($cmd);
Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;


##
print "get potential truncation position ...\n" if $verbose;
my $verboseFlag = $verbose ? '-v' : '';
my $uniqTagTruncBedFile = "$cache/tag.trunc.bed";

$cmd = "perl $cmdDir/bedExt.pl $verboseFlag -n up -l \"-1\" -r \"-1\" $uniqTagCleanBedFile $uniqTagTruncBedFile";
print $cmd, "\n" if $verbose;

$ret = system($cmd);
Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;

##
print "cluster unique tags...\n" if $verbose;
my $bigFlag = $bigFile ? '-big' : '';
my $keepCacheFlag = $keepCache ? '--keep-cache' : '';
my $tagClusterBedFile = "$cache/tag.cluster.0.bed";

$cmd = "perl $cmdDir/tag2cluster.pl $bigFlag $verboseFlag $keepCacheFlag -c $cache/tag2cluster_cache -s -maxgap \"-1\" $uniqTagBedFile $tagClusterBedFile";
print $cmd, "\n" if $verbose;

$ret = system($cmd);
Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;

my $tagClusterCleanBedFile = "$cache/tag.cluster.bed";
$cmd = "awk '{if(\$5>2) {print \$0}}' $tagClusterBedFile > $tagClusterCleanBedFile";
print $cmd, "\n" if $verbose;

$ret = system($cmd);
Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;

##

print "identifying CITS ...\n" if $verbose;
my $multTestFlag = $multiTestCorrection ? '--multi-test' : '';

$cmd = "perl $cmdDir/tag2peak.pl $bigFlag $verboseFlag $keepCacheFlag -c $cache/tag2peak_cache -ss --prefix $prefix -gap $maxGap -p $pvalueThreshold $multTestFlag --gene $tagClusterCleanBedFile $uniqTagTruncBedFile $outCITSBedFile";
print $cmd, "\n" if $verbose;

$ret = system($cmd);
Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;

system ("rm -rf $cache") unless $keepCache;


