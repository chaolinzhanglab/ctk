#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Carp;
use MyConfig;

my $prog = basename ($0);
my $cmdDir = dirname ($0);

my $verbose = 0;
my $cache = getDefaultCache ($prog);


GetOptions ("cache:s"=>\$cache,
		"v"=>\$verbose);


if (@ARGV != 3)
{
	print "get CIMS allele frequency\n";
	print "Usage: $prog [options] <mutation.txt> <CIMS.bed> <out.txt>\n";
	print " -cache: cache dir ($cache)\n";
	print " -v    : verbose\n";
	exit (1);
}


my ($inMutationFile, $inCIMSFile, $outFile) = @ARGV;


Carp::croak "$cache already exists\n" if -d $cache;
system("mkdir $cache");

my %CIMStags;
my %CIMSHash; 
my $cmd;

# 1. get the number of tags at each mutation position
print "Counting mutations ...\n" if $verbose;
open (INFILE, $inCIMSFile) || Carp::croak "can not open file $inCIMSFile to read\n";
while (my $line = <INFILE>)
{
	chomp $line;
	my @cols = split ("\t", $line);
	my $clusterId = $cols[3];
	$clusterId =~ m/\[k=(\d+)\]/;
	$CIMStags{$clusterId} = $1;
}

close(INFILE);


# 2. get allele frequency from mutation file
# 2a. reformat mutation file into bed (and reverse complement)
print "Process mutations ...\n" if $verbose;
my $formatMutationFile = "$cache/mutation.sub.bed";
open (OUTFILE, ">$formatMutationFile") || Carp::croak "can not open file $formatMutationFile to write\n";
open (INFILE, "<$inMutationFile") || Carp::croak "can not open file $inMutationFile to read\n";
while (my $line = <INFILE>)
{
	my @cols = split ("\t", $line);
	next if $cols[8] ne ">"; # check that it is a substitution
	my $refBase = $cols[7];
	my $subBase = $cols[9];

	if ($cols[5] eq "-") # reverse complement minus strand
	{
		$refBase =~ tr/ATGC/TACG/;
		$subBase =~ tr/ATGC/TACG/;
	}

	print OUTFILE "$cols[0]\t$cols[1]\t$cols[2]\t$cols[3]^^^$refBase^$subBase\t$cols[4]\t$cols[5]\n"; # print with ref/mutation bases added to tag name
}
close (INFILE);
close (OUTFILE);

# 3. overlap with CIMS
print "compare mutations and CIMS ...\n" if $verbose;
my $mutationVsCIMS = "$cache/mutation.vs.CIMS";
$cmd = "perl $cmdDir/tagoverlap.pl -region $inCIMSFile -ss $formatMutationFile $mutationVsCIMS";
$cmd .= " -v" if $verbose;
system ($cmd);


# 4. count alleles
# 4a. retrieve mutations
print "parse tag/cluster comparison ...\n" if $verbose;

open (INFILE, "<$mutationVsCIMS") || Carp::croak "can not open file $mutationVsCIMS to read\n";



while (my $line = <INFILE>)
{
	next if $line =~/^\s*$/;

	my @cols = split ("\t", $line);
	my $name = $cols[3];
	my $strand = $cols[5];

	my ($tagId, $clusterId) = split (/\/\//, $name);
	
	#$clusterId =~/\/(\d+)\]/;
	#my $total = $1;

	#$CIMSHash{$clusterId}->{'total'} = $total;
	my $refAllele = substr ($tagId, -3, 1);
	my $subAllele = substr ($tagId, -1, 1);
	$subAllele = uc ($subAllele);
	#$subAllele =~tr/ACGT/TGCA/ if $strand eq '-';
	$CIMSHash{$clusterId}->{$subAllele}++;
	$CIMSHash{$clusterId}->{"ref"} = uc($refAllele);
}

close (INFILE);

# 4b. count mutations
print "Count mutations ...\n" if $verbose;

my $outAlleleFile = "$cache/allele.txt";
open (OUTFILE, ">$outAlleleFile") || Carp::croak "cannot open file $outFile to write\n";
print OUTFILE join ("\t", "CIMS", "A", "C", "G", "T"), "\n";
foreach my $clusterId (sort keys %CIMSHash)
{
	my $cluster = $CIMSHash{$clusterId};
	my $refAllele = $cluster->{"ref"};

	if (!exists $cluster->{'A'}) {$cluster->{'A'} = 0};
	if (!exists $cluster->{'C'}) {$cluster->{'C'} = 0};
	if (!exists $cluster->{'G'}) {$cluster->{'G'} = 0};
	if (!exists $cluster->{'T'}) {$cluster->{'T'} = 0};

	
	if ($cluster->{$refAllele} != 0 ) 
	{
		Carp::croak "for cluster $cluster, counts for reference allele $refAllele are not 0\n";
	}
	else
	{
		$cluster->{$refAllele} = $CIMStags{$clusterId} - $cluster->{'A'} - $cluster->{'T'} - $cluster->{'G'} - $cluster->{'C'};
	}
	
	my $A = $cluster->{'A'};
	my $C = $cluster->{'C'};
	my $G = $cluster->{'G'};
	my $T = $cluster->{'T'};

	print OUTFILE join ("\t", $clusterId, $A, $C, $G, $T), "\n";
}
close (OUTFILE);

$cmd = "perl $cmdDir/selectRow.pl -f 3 $outAlleleFile $inCIMSFile > $outFile";
system ($cmd);

print "Done.\n";
