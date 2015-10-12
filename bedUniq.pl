#!/usr/bin/perl -w

use lib "/data/zhangc/perl_lib2";

use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use Carp;

use Bed;

my $prog = basename($0);

#TO DO: replace subroutines with those in Bed.pm

my $sameStrand = 0;	#cluster regions are in the same strand for a match
my $sameName = 0; #cluster only regions with the same name
my $maxGap = 0; #the max gap to be considered as an overlap
my $criterion = 'maxScore'; #or 'minScore', 'maxLen', 'minLen', 'random', 'singleton'
#singleton means no overlap with other regions

my $strictOverlap = 0; #exactly the same interval
my $debug = 0;
my $verbose = 0;

GetOptions ('c|criteria:s'=>\$criterion,
			'so|same-interval'=>\$strictOverlap,
			'ss|same-strand'=>\$sameStrand,
			'sn|same-name'=>\$sameName,
			'maxgap:i'=>\$maxGap,
			'd|debug'=>\$debug,
			'v|verbose'=>\$verbose
			);


if (@ARGV != 2)
{
	print "get nonoverlapping unique regions\n";
	
	print "Usage: $prog [options] <in.bed> <out.bed>\n";
	print " <in.bed> and <out.bed>: use - for STDIN and STDOUT\n";
	print " -c: criterion ([maxScore]|minScore|maxLen|minLen|random|singleton)\n";
	print " -so: same interval required [off]\n";
	print " -ss: same strand required [off]\n";
	print " -sn: same name required [off]\n";
	print " -maxgap: the max gap to be considered as an overlap\n";
	print " -d: debug (on|[off])\n";
	print " -v: verbose (on|[off])\n";
	exit (1);
}

print "region selection criterion: $criterion\n" if $verbose;
Carp::croak "maxgap can not be negative: $maxGap\n" if $maxGap < 0;


my ($inBedFile, $outBedFile) = @ARGV;

print "loading regions from $inBedFile ...\n" if $verbose;
my $inRegions = readBedFile ($inBedFile, $verbose);
my $n = @$inRegions;
print "$n regions loaded\n" if $verbose;


Carp::croak "No strand specified\n" if $sameStrand && (not exists $inRegions->[0]->{"strand"});
Carp::croak "No name specified\n" if $sameName && (not exists $inRegions->[0]->{"name"});


my $regionIdx;
foreach my $r (@$inRegions)
{
	#my $chrom = $r->{"chrom"};
	
	my $key = $r->{"chrom"};
	$key .= "//" . $r->{"strand"} if ($sameStrand);
	$key .= "//" . $r->{"name"} if ($sameName);
	
	push @{$regionIdx->{$key}}, $r;

	my $n = @{$regionIdx->{$key}};
	$regionIdx->{$key}->[$n - 1]->{"idx"} = $n - 1;
}
my $nkey = keys %$regionIdx;
print "clustering $nkey groups separately\n" if $verbose;

sub clusterRegionsExact
{
	my $regionsOnChrom = $_[0];
	my %regionHash;
	foreach my $r (@$regionsOnChrom)
	{
		my $key = $r->{"chrom"} . ":" . $r->{"chromStart"} . "-" . $r->{"chromEnd"};
		push @{$regionHash{$key}}, $r->{"idx"};
	}
	my @clusters;
	foreach my $key (keys %regionHash)
	{
		push @clusters, $regionHash{$key};
	}
	return \@clusters;
}


sub clusterRegions2
{
	my $regionsOnChrom = $_[0]; 
	#only regions in the same group. they are on the same chromosome, 
	#but be filtered by strand and/or name
	
	#sort ascendingly according to the start point
	my @regionsOnChromSorted = sort {$a->{"chromStart"} <=> $b->{"chromStart"}} @$regionsOnChrom;
	
	my $n = @regionsOnChromSorted;

	print "\n\nclustering begins with $n regions...\n" if $debug;
	
	my @clusters;
	
	my @currCluster;
	my $currClusterEnd = -$maxGap - 1;
	my $i = 0;
	foreach my $r (@regionsOnChromSorted)
	{
		#next if ($strand ne 'b' && $strand ne $r->{"strand"});
		$i++;

		my $chromStart = $r->{"chromStart"} - $maxGap / 2;
		my $chromEnd = $r->{"chromEnd"} + $maxGap / 2;
		my $idx = $r->{"idx"};
		
		print "$i: clusterEnd = $currClusterEnd, chromStart = $chromStart, chromEnd = $chromEnd, idx=$idx\n" if $debug;
		#begin a new cluster
		if ($chromStart > $currClusterEnd)
		{
			my $n = @currCluster;
			print "close the old cluster with $n regions ...\n" if $debug;
			print join ("\t", @currCluster), "\n" if $debug;

			if (@currCluster > 0)
			{
				my @currClusterCpy = @currCluster;
				push @clusters, \@currClusterCpy;
			}
			print "begin a new cluster...\n" if $debug;
			@currCluster = ();
			push @currCluster, $r->{"idx"};
			$currClusterEnd = $chromEnd;
		}
		else
		{
			my $n = @currCluster;
			print "expand the current cluster with $n regions ...\n" if $debug;
			#expand the old cluster
			push @currCluster, $r->{"idx"};
			$currClusterEnd = $chromEnd if $currClusterEnd < $chromEnd;
		}
	}
	
	if (@currCluster >0)
	{
		my @currClusterCpy = @currCluster;
		push @clusters, \@currClusterCpy;
	}
	
	$n = @clusters;
	print "\n\n $n clusters found\n\n" if $debug;
	return \@clusters;
}


sub selectRegions
{
	my ($regions, $clusters) = @_;
	
	my $nc = @$clusters;
	
	foreach my $c (@$clusters)
	{
		my $n = @$c;
		print "\ncluster $c: $n regions\n" if $debug;
		

		my $str = "";
		foreach my $idx (@$c)
		{
			my $r = $regions->[$idx];
			$str .= "//" . $r->{"name"};
		}
		
		if ($debug)
		{
			foreach my $idx (@$c)
			{
				#print "idx=$idx\n";
				my $r = $regions->[$idx];
				my $score = exists $r->{"score"} ? $r->{"score"} : 0;
				print join ("\t", $r->{"chrom"}, $r->{"chromStart"}, $r->{"chromEnd"}, $score), "\n";
			}
		}
		if ($n == 1)
		{
			my $idx = $c->[0];
			my $r = $regions->[$idx];
			$r->{'print'} = 1;
			next;
		}

		if ($criterion eq 'singleton')
		{
			next;
		}
		
		if ($criterion eq 'random')
		{
			my $i = int(rand ($n));
			my $idx = $c->[$i];
			my $r = $regions->[$idx];
			$r->{'print'} = 1;

			if ($debug)
			{
				print "-->selected region:\n";
				my $score = exists $r->{"score"} ? $r->{"score"} : 0;
				print join ("\t", $r->{"chrom"}, $r->{"chromStart"}, $r->{"chromEnd"}, $score), "\n";
			}
			next;
		}

		my @regionsInCluster;
		my @regionSorted;
		foreach my $idx (@$c)
		{
			my $r = $regions->[$idx];
			$r->{'len'} = $r->{'chromEnd'} - $r->{'chromStart'} + 1;
			push @regionsInCluster, $r;
		}
		
		if ($criterion eq 'maxScore')
		{
			@regionSorted = sort {$b->{'score'} <=> $a->{'score'}} @regionsInCluster;
		}
		elsif ($criterion eq 'minScore')
		{
			@regionSorted = sort {$a->{'score'} <=> $b->{'score'}} @regionsInCluster;
		}
		elsif ($criterion eq 'maxLen')
		{
			@regionSorted = sort {$b->{'len'} <=> $a->{'len'}} @regionsInCluster;
		}
		elsif ($criterion eq 'minLen')
		{
			@regionSorted = sort {$a->{'len'} <=> $b->{'len'}} @regionsInCluster;
		}
		else
		{
			Carp:croak "invalid value for selection criterion: $criterion\n";
		}
		
		
		my $r = $regionSorted[0];
		$r->{"print"} = 1;
		#$r->{"name"} = $str;
		if ($debug)
		{
			print "-->selected regions:\n";
			my $score = exists $r->{"score"} ? $r->{"score"} : 0;
			print join ("\t", $r->{"chrom"}, $r->{"chromStart"}, $r->{"chromEnd"}, $score), "\n";
		}
	}
}



foreach my $key (sort keys %$regionIdx)
{
	my $regionsInGroup = $regionIdx->{$key};
	
	print "\n clustering regions in group $key...\n" if $debug;
	my $clusters;
   	if ($strictOverlap)
	{
		$clusters = clusterRegionsExact ($regionsInGroup);
	}
	else
	{
		$clusters = clusterRegions2 ($regionsInGroup);
	}	
	my $nc = @$clusters;
	print "\n\n$nc clusters found\n" if $debug;
	selectRegions ($regionsInGroup, $clusters);
}

my @outRegions;
foreach my $r (@$inRegions)
{
	push @outRegions, $r if exists $r->{"print"};
}

writeBedFile (\@outRegions, $outBedFile);





