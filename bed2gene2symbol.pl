#!/usr/bin/perl -w

use strict;

use File::Basename;
use MyConfig;
use Getopt::Long;

my $prog = basename ($0);
my $cmdDir = dirname ($0);

my $cache = getDefaultCache($prog);
my $verbose = 0;
my $separateStrand = 0;
my $noRegionId = 0;
my $big = 0;

GetOptions ("cache:s"=>\$cache,
		"ss"=>\$separateStrand,
		"no-region-id"=>\$noRegionId,
		"big"=>\$big,
		"v"=>\$verbose);



if (@ARGV != 5)
{
	print "mapping regions to genes and symbols\n";
	print "Usage: $prog [options] <in.bed> <refGene_knownGene.bed> <refGene_knownGene_To_gene.txt> <gene2symbol.txt> <out.txt>\n";
	print " -ss             : separate strand\n";
	print " -big            : big file\n";
	print " --no-region-id  : do not print region id\n";
	print " -cache  [string]: cache dir ($cache)\n";
	print " -v              : verbose\n";
	exit (0);
}


my ($inBedFile, $transcriptBedFile, $transcriptToGeneFile, $gene2SymbolFile, $outFile) = @ARGV;


system ("mkdir $cache");

my $bigFlag = $big ? '-big' : '';
my $verboseFlag = $verbose ? '-v' : '';

my $tagVSTs = "$cache/in.vs.ts.bed";
my $ssFlag = $separateStrand ? "-ss" : " ";
system ("perl $cmdDir/tagoverlap.pl $bigFlag $ssFlag $verboseFlag -region $transcriptBedFile $inBedFile $tagVSTs");
Carp::croak "cmd failed\n" unless -f $tagVSTs;

my $tagVSTsIdpair = "$cache/in.vs.ts.idpair";
my $cmd = "awk \'{print \$4}\' $tagVSTs | awk -F \"//\" \'{print \$1\"\\t\"\$2}\' > $tagVSTsIdpair";
system ($cmd);
Carp::croak "cmd=$cmd failed\n" unless -f $tagVSTsIdpair;

my $tagVSTs_Tsid = "$cache/in.vs.ts.tsid";
$cmd = "awk \'{print \$2}\' $tagVSTsIdpair > $tagVSTs_Tsid";
system ($cmd);

my $tsidToGeneId = "$cache/tsid2geneid";
$cmd = "perl $cmdDir/selectRow.pl -p -pt \"\" -s $transcriptToGeneFile $tagVSTs_Tsid > $tsidToGeneId";
system ($cmd);
Carp::croak "cmd=$cmd failed\n" unless -f $tsidToGeneId;

my $tag2gene = "$cache/tag2gene";
$cmd = "paste $tagVSTsIdpair $tsidToGeneId | awk \'{if (NF>3) {print \$1\"\\t\"\$4}}' | sort | uniq > $tag2gene";
system ($cmd);
Carp::croak "cmd=$cmd failed\n" unless -f $tag2gene;

my $tag2geneUniq = "$cache/tag2gene.uniq";
$cmd = "perl $cmdDir/uniqRow.pl $verboseFlag $tag2gene $tag2geneUniq";
system ($cmd);
Carp::croak "cmd=$cmd failed\n" unless -f $tag2geneUniq;

my $tag2geneUniq_geneid = "$cache/tag2gene.uniq.geneid";
$cmd = "awk \'{print \$2}\' $tag2geneUniq > $tag2geneUniq_geneid";
system ($cmd);

my $gene2symbol = "$cache/gene2symbol";
$cmd = "perl $cmdDir/selectRow.pl -p -pt \"\" -s $gene2SymbolFile $tag2geneUniq_geneid > $gene2symbol";
system ($cmd);
Carp::croak "cmd=$cmd failed\n" unless -f $gene2symbol;


my $tag2gene2symbol = "$cache/tag2gene2symbol";
$cmd = "paste $tag2geneUniq $gene2symbol | awk \'{print \$1\"\\t\"\$3\"\\t\"\$4}\' > $tag2gene2symbol";
system ($cmd);
Carp::croak "cmd=$cmd failed\n" unless -f $tag2gene2symbol;

my $tagid = "$cache/tagid";

$cmd = "awk \'{print \$4}\' $inBedFile > $tagid";
system ($cmd);
Carp::croak "cmd=$cmd failed\n" unless -f $tagid;


$cmd = "perl $cmdDir/selectRow.pl -p -pt \"\" -s $tag2gene2symbol $tagid | awk '{print \$1\"\\t\"\$2\"\\t\"\$3}' > $outFile";
$cmd = "perl $cmdDir/selectRow.pl -p -pt \"\" -s $tag2gene2symbol $tagid | awk '{print \$2\"\\t\"\$3}' > $outFile" if $noRegionId;

system ($cmd);

system ("rm -rf $cache");

