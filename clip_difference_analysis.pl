#!/usr/bin/perl -w

use strict;
use warnings;
use Carp;
use Getopt::Long;
use File::Basename;
use Math::CDF qw(:all);


use Bed;
use Common;

my $verbose = 0;
my $prog = basename ($0);
my $filterIDFile = "";
my $pseudoCount = 1e-7;

my $exprFile = "";


GetOptions ("v|verbose"=>\$verbose,
	'filter:s'=>\$filterIDFile,
	'pseudo-count:f'=>\$pseudoCount,
	'expr:s'=>\$exprFile,
	);

if (@ARGV != 3)
{
	print "CLIP cluster difference analysis\n";
	print "Usage: $prog <clip1.bed> <clip2.bed> <out.txt>\n";
	print " -filter       [string] : list of clusters to estiamte expected ratio\n";
	print " -pseudo-count [float]  : pseudo count of tag number ($pseudoCount)\n";
	print " -expr         [string] : expression file (cluster id<\\t>log2fc\n";
	print " -v                     : verbose\n";
	exit (1);
}


my ($clusterBedFile1, $clusterBedFile2, $outFile) = @ARGV;

print "read clusters from $clusterBedFile1 ...\n" if $verbose;
my $clustersCondition1 = readBedFile ($clusterBedFile1, $verbose); 

print "read clusters from $clusterBedFile2 ...\n" if $verbose;
my $clustersCondition2 = readBedFile ($clusterBedFile2, $verbose); 

my $n = @$clustersCondition1;

Carp::croak "the two input BED files have different number of rows\n" unless $n == @$clustersCondition2;
print "$n clusters loaded\n" if $verbose;


my %filterIDHash;
if ($filterIDFile ne '')
{
	Carp::croak "$filterIDFile does not exist\n" unless -f $filterIDFile;

	print "reading cluster ids to estimate expected ratio...\n" if $verbose;
	my $fin;
	open ($fin, "<$filterIDFile") || Carp::croak "cannot open file $filterIDFile ...\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;
		my @cols = split(/\s/, $line);

		$filterIDHash{$cols[0]} = 1;
	}
	close ($fin);

	my $n = keys %filterIDHash;
	print "$n clusters loaded\n" if $verbose;
}

my $totalTagNumCondition1 = 0;
my $totalTagNumCondition2 = 0;

for (my $i = 0; $i < $n; $i++)
{
	my $name = $clustersCondition1->[$i]->{'name'};
	Carp::croak "cluster $i in the two condition does not have the same name:\n" , Dumper ($clustersCondition1->[$i]), Dumper ($clustersCondition2->[$i]), "\n"
	unless $clustersCondition2->[$i]->{'name'} eq $name;

	if ($filterIDFile ne '')
	{
		next unless exists $filterIDHash{$name};
	}

	$totalTagNumCondition1 += $clustersCondition1->[$i]->{'score'};
	$totalTagNumCondition2 += $clustersCondition2->[$i]->{'score'};
}


my $expectedRatio = $totalTagNumCondition1 / $totalTagNumCondition2;

print "total tag number N1=$totalTagNumCondition1, N2=$totalTagNumCondition2, N1/N2=$expectedRatio\n" if $verbose;

my @exprLog2FC;

if ($exprFile ne '' && (-f $exprFile))
{
	print "reading gene expression fold change ...\n" if $verbose;

	my $fin;
	open ($fin, "<$exprFile") || Carp::croak "cannot open file $exprFile to read\n";
	my $iter = 0;

	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line=~/^\s*$/;

		print "$iter ...\n" if $verbose && $iter % 100000 == 0;
		my ($clusterId, $log2FC) = split (/\s/, $line);
		
		Carp::croak "inconsistency detected in cluster ids of input files\n" unless $clusterId eq $clustersCondition1->[$iter]->{'name'} && $clusterId eq $clustersCondition2->[$iter]->{'name'};
		push @exprLog2FC, $log2FC;

		$iter++;
	}
	close ($fin);
}

if (@exprLog2FC > 0)
{
	Carp::croak "expression file and cluster bed file must have the same number of lines in the same order\n" unless @exprLog2FC == $n;
}

print "performing binomial test of difference\n" if $verbose;

#&R::initR("--silent");
#&R::library("RSPerl");

my $fout;
open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile\n";
print $fout join ("\t", "chrom", "chromStart", "chromEnd", "name", "score", "strand", "n1", "n2", "log2FC", "p"), "\n";
for (my $i = 0; $i < $n; $i++)
{
	print "$i ...\n" if $verbose && $i % 100000 == 0;

	my $n1 = $clustersCondition1->[$i]->{'score'} > 0 ? $clustersCondition1->[$i]->{'score'} : $clustersCondition1->[$i]->{'score'} + $pseudoCount; #1 is the pseudo count
	my $n2 = $clustersCondition2->[$i]->{'score'} > 0 ? $clustersCondition2->[$i]->{'score'} : $clustersCondition2->[$i]->{'score'} + $pseudoCount;

	my $r = $expectedRatio;
	if (@exprLog2FC == $n)
	{
		my $log2FC = $exprLog2FC[$i];
		my $fc = 2**$log2FC;
		$r *= $fc;
	}

	my $p = $r / ($r + 1);

	my $pvalue = binomTest (int($n1+0.5), int($n1+$n2+0.5), $p);
	#my $ret = &R::call ("binom.test", int($n1+0.5)+0, int($n1+$n2+0.5)+0, $p + 0);
	#my $pvalue = $ret->getEl ('p.value');

	#$pvalue = 1- $pvalue if $pvalue > 0.5;
	#$pvalue *= 2; #make it two tails

	$pvalue = 1e-100 if $pvalue <= 0;
	$pvalue = 1 if $pvalue > 1;

	my $log2FC = log ($n1 / $n2 / $r) / log(2);
	
	my %c = %{$clustersCondition1->[$i]};
	$c{'score'} = sprintf ("%.2f", -log ($pvalue) / log(10) * 10); 

	print $fout join ("\t", bedToLine (\%c), $n1, $n2, $log2FC, $pvalue), "\n";
}

close ($fout);




