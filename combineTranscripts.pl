#!/usr/bin/perl


use strict;
use warnings;

use Getopt::Long;
use File::Basename;

use Bed;


my $prog = basename ($0);
my $big = 0;

my $verbose = 0;
my $inconsistBedFile = "";
my @ARGV0 = @ARGV;

GetOptions ("print-inconsist:s"=>\$inconsistBedFile,
		"v|verbose"=>\$verbose);

if (@ARGV != 3)
{
	print "combine transcripts of the same gene to get the exonic region\n";
	print "$prog [options] <transcript.bed> <ts2gene> <gene.exonic.bed>\n";
	print " --print-inconsist  [string]: print transcripts with inconsistency\n";
	print " -v               : verbose\n";
	exit (1);
}


my ($transcriptBedFile, $transcript2geneFile, $combinedBedFile) = @ARGV;

print "CMD=$prog ", join(" ", @ARGV0), "\n" if $verbose;

print "reading transcripts from $transcriptBedFile ...\n" if $verbose;

my $transcripts = readBedFile ($transcriptBedFile, $verbose);
my $n = @$transcripts;

print "$n transcripts loaded\n" if $verbose;

my %transcriptHash;

foreach my $r (@$transcripts)
{
	my $tsId = $r->{"name"};
	#$transcriptHash{$tsId} = $r;
	push @{$transcriptHash{$tsId}}, $r; #in case a transcript mapped to multiple loci or represent exons
}

print "reading transcript to gene mapping from $transcript2geneFile ...\n" if $verbose;


my $fin;
my %geneHash;
open ($fin, "<$transcript2geneFile") || Carp::croak "can not open file $transcript2geneFile to read\n";
while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;
	my ($tsId, $geneId) = split(/\t/, $line);
	push @{$geneHash {$geneId}}, $tsId;
}
close ($fin);

$n = keys %geneHash;

print "$n genes loaded\n" if $verbose;


print "combine transcripts for each gene...\n" if $verbose;

my $iter = 0;
my (@genes, @remainingTranscripts);
foreach my $geneId (sort keys %geneHash)
{
	$iter++;

	my $transcripts_gene = $geneHash{$geneId};
	my @regions;

	foreach my $tsId (@$transcripts_gene)
	{
		next unless exists $transcriptHash{$tsId};
		#push @regions, $transcriptHash{$tsId};
		push @regions, @{$transcriptHash{$tsId}};
	}

	next unless @regions > 0;

	my $chrom = $regions[0]->{"chrom"};
	my $strand = $regions[0]->{"strand"};

	my $consist = 1;
	for (my $i = 1; $i < @regions; $i++)
	{
		my $r = $regions[$i];
		if ($chrom ne $r->{"chrom"} || $strand ne $r->{"strand"})
		{
			$consist = 0;
			print "inconsistency in gene $geneId found\n";
			last;
		}
	}

	if ($consist)
	{
		my $g = combineRegions (\@regions);
		$g->{"name"} = $geneId;
		push @genes, $g;
	}
	else
	{
		push @remainingTranscripts, @regions;
	}
}

print "dump results to $combinedBedFile File ...\n" if $verbose;

writeBedFile (\@genes, $combinedBedFile);
writeBedFile (\@remainingTranscripts, $inconsistBedFile) if $inconsistBedFile ne '';


=end
if (0)
{
	foreach my $gene (@genes)
	{
		my $length = Common::sum ($gene->{"blockSizes"});
		print join ("\t", $gene->{"name"}, $length), "\n";
	}
}

