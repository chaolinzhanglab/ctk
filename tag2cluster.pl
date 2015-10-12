#!/usr/bin/perl -w

use strict;
use warnings;

use File::Basename;
use Data::Dumper;
use Getopt::Long;
use Carp;

use MyConfig;
use Bed;

my $prog = basename($0);

my $bigFile = 0;	#if yes, we need to use cache

my $sameStrand = 0;	#cluster regions are in the same strand for a match
my $maxGap = 0; #the max gap to be considered as an overlap
my $overlapFraction = 0; #1e-9; #if we don't allow gap, what the fraction of overlap we consider a match
my $collapseMode = 0; #0-- no collapse; 1 -- collapse if both ends match; 2-- collapse if one falls in another
my $randomLinker = 0; #do not collapse if the linker sequences are different
my $keepScore = 0; #keep the original score of representative tag (for collapse mode)

#singleton means no overlap with other regions
my $weight = 0;	#the way to calculate the score in each cluster, without weight means the number of tags, other wise, the sum of scores
my $weightInName = 0;

my $debug = 0;
my $verbose = 0;
my $outputFormat = 'bed'; #or wig

my $cache = getDefaultCache ($prog);
my $keepCache = 0;

my @ARGV0 = @ARGV;

GetOptions (
			'big'=>\$bigFile,
			's|same-strand'=>\$sameStrand,
			'weight'=>\$weight,
			'weight-in-name'=>\$weightInName,
			'maxgap:i'=>\$maxGap,
			'collapse:i'=>\$collapseMode,
			'overlap:f'=>\$overlapFraction,
			'rl|random-linker'=>\$randomLinker,
			'keep-score'=>\$keepScore,
			'of:s'=>\$outputFormat,
			'c|cache:s'=>\$cache,
			'keep-cache'=>\$keepCache,
			'd|debug'=>\$debug,
			'v|verbose'=>\$verbose
			);


if (@ARGV != 2)
{
	print "cluster overlapping tags into clusters and report tag count in each cluster\n";
	
	print "Usage: $prog [options] <in.bed> <out.bed>\n";
	print " <in.bed>         : gzip compressed file with .gz extension is acceptable\n";
	print "                    use \"-\" for stdin\n";
	print "options\n";
	print " -big             : set when the input file is big\n";
	print " -s               : same strand required [off]\n";
	print " -weight          : consider the weight of each tag\n";
	print " --weight-in-name : find weight in name\n";
	print " -maxgap    [int] : the max gap to be considered as an overlap ($maxGap)\n";
	print " -collapse  [int] : collapse mode (0: no collapse; 1: collapse if match both ends; 2: collapse if one is in another)\n";
	print " -overlap [float] : overlap fraction, effective unless collapse ($overlapFraction)\n";
	#print " --random-linker  : random linker exists, no collapse for different linkers\n";
	#print " --keep-score     : keep the original score of representative tags (for collapse mode)\n";
	print " -of     [string] : output format ([bed]|wig)\n";
	print " -c      [string] : cache dir ($cache)\n";
	print " --keep-cache     : keep cache when the job done\n";
	print " -d               : debug (on|[off])\n";
	print " -v               : verbose (on|[off])\n";
	exit (1);
}

print "CMD=$prog ", join(" ", @ARGV0), "\n" if $verbose;

#Carp::croak "maxgap can not be negative: $maxGap\n" if $maxGap < 0;
Carp::croak "collapse must be 0, 1 or 2\n" if $collapseMode * ($collapseMode -1) * ($collapseMode - 2) != 0;

$maxGap = 0 if $collapseMode;


my ($inBedFile, $outBedFile) = @ARGV;

if ($collapseMode)
{
	Carp::croak "output format must be bed\n" if $outputFormat ne 'bed';
}



Carp::croak "$cache already exists\n" if -d $cache;
system ("mkdir $cache");


my %tagCount;

print "reading tags from $inBedFile ...\n" if $verbose;
if ($bigFile)
{
	my $ret = splitBedFileByChrom ($inBedFile, $cache, $verbose);
	%tagCount = %$ret;
}
else
{
	my $tags = readBedFile ($inBedFile, $verbose);
	foreach my $t (@$tags)
	{
		my $chrom = $t->{"chrom"};
		push @{$tagCount{$chrom}}, $t;
	}
}

print "get tag count broken down into chromosomes ...\n" if $verbose;

foreach my $chrom (sort keys %tagCount)
{

	my $n = $tagCount{$chrom};
	$n = ref($n) eq 'HASH' ? $n = $n->{'n'} : @$n;

	print "$chrom : $n tags\n" if $verbose;
}


print "\n\nclustering tags ...\n" if $verbose;

my @strands = ('b');
@strands = qw(+ -) if $sameStrand;


my $fout;
if ($outBedFile eq "-")
{
     $fout = *STDOUT;
}
else
{
	open ($fout, ">$outBedFile") || Carp::croak "can not open file $outBedFile to write\n";
}

foreach my $s (@strands)
{
	
	print "processing strand $s ...\n" if $verbose;
	
	if ($outputFormat eq 'wig')
	{
		print $fout "\n\ntrack type=wiggle_0 name=tag2cluster$s autoscale=on\n";
	}

	foreach my $chrom (sort keys %tagCount)
	{
		my $tags;
		if ($bigFile)
		{
			my $tmpFile = $tagCount{$chrom}->{'f'};
			print "loading tags on chromsome $chrom from $tmpFile...\n" if $verbose;
			$tags = readBedFile ($tmpFile, $verbose);
		}
		else
		{
			$tags = $tagCount{$chrom};
		}

		my $n = @$tags;

		if ($weight && $weightInName)
		{
			foreach my $t (@$tags)
			{
				my @cols = split (/\#/, $t->{'name'});
				pop @cols if $randomLinker;
				$t->{'score'} = pop @cols;
			}
		}

		print "$n tags loaded on chromosome $chrom\n" if $verbose;
	
		Carp::croak "No strand specified\n" if $sameStrand && (not exists $tags->[0]->{"strand"});

		print "clustering $chrom ...\n" if $verbose;
		#my @tags4cluster;
		#my $i = 0;
		#foreach my $t (@$tags)
		#{
		#	next if $s ne 'b' && $t->{"strand"} ne $s;
		#	$t->{"idx"} = $i++;
		#	push @tags4cluster, $t;
		#}
		#
		#my $clusters = clusterRegions (\@tags4cluster, $s);

		my $clusters;
	   
		if ($collapseMode && $randomLinker)
		{
			#because we need to run clustering again in cluster mode, using only tags with the same linker
			#definition clusterRegions ($regionsOnChrom, $strand, $maxGap, $overlapFraction, $collapse); 
			
			$clusters = clusterRegions ($tags, $s, 0, 0, 0);
		}
		else
		{
			$clusters = clusterRegions ($tags, $s, $maxGap, $overlapFraction, $collapseMode);
		}
		my $nc = @$clusters;
		print "\n\n$nc clusters found\n" if $debug;

		my $iter = 0;
		foreach my $clust (@$clusters)
		{
			#my @tagsInClust = @tags4cluster[@$clust];
			my @tagsInClust = @$tags[@$clust];
			my $score = $weight ? 0 : @tagsInClust;

			my $chromStart = 1e9;
			my $chromEnd = -1;
			foreach my $t (@tagsInClust)
			{
				if ($weight)
				{
					Carp::croak "not score found in tag:", Dumper ($t), "\n" unless exists $t->{"score"};
					$score += $t->{"score"};
				}
				$chromStart = $t->{"chromStart"} if $t->{"chromStart"} < $chromStart;
				$chromEnd = $t->{"chromEnd"} if $t->{"chromEnd"} > $chromEnd;
			}

			#my $score = @tagsInClust;

			
			my $name = $chrom;
		    $name .= "_f_c$iter" if $s eq '+';
			$name .= "_r_c$iter" if $s eq '-';
			
			$name .= "_c$iter" if $s eq 'b';
		
			if ($outputFormat eq 'wig')
			{
				print $fout join ("\t", $chrom, $chromStart, $chromEnd + 1, $score), "\n";
				next;
			}
			
			if ($collapseMode)
			{
				@tagsInClust = sort {($b->{"chromEnd"} - $b->{"chromStart"}) <=> ($a->{"chromEnd"} - $a->{"chromStart"})} @tagsInClust;

				if ($randomLinker)
				{
					my %tagsInClust;
					#divide according to linker sequence
					foreach my $tag (@tagsInClust)
					{
						my $tagName = $tag->{"name"};
						my @cols = split (/\#/, $tagName);
						my $linker = pop @cols;
						next if $linker =~/[^ACGT]/i;

						push @{$tagsInClust{$linker}}, $tag;
					}


					foreach my $linker (keys %tagsInClust)
					{
						my $tagsInClustWithSameLinker = $tagsInClust{$linker};

						#cluster again for tags with the same linker in the cluster
						my $clusters2 = clusterRegions ($tagsInClustWithSameLinker, $s, $maxGap, $overlapFraction, $collapseMode);

						foreach my $clust2 (@$clusters2)
						{
							my @tagsInClust2 = @$tagsInClustWithSameLinker[@$clust2];
							my $score = $weight ? 0 : @tagsInClust2;

							foreach my $t (@tagsInClust2)
							{
								if ($weight)
								{
									Carp::croak "not score found in tag:", Dumper ($t), "\n" unless exists $t->{"score"};
									$score += $t->{"score"};
								}
							}
							$tagsInClust2[0]->{"score"} = $score unless $keepScore;
							print $fout bedToLine ($tagsInClust2[0]), "\n";
						}
					}
				}
				else
				{
					$tagsInClust[0]->{"score"} = $score unless $keepScore;
					print $fout bedToLine ($tagsInClust[0]), "\n";
				}
			}	
			else
			{
				if ($s ne 'b')
				{
					print $fout join ("\t", $chrom, $chromStart, $chromEnd + 1, $name, $score, $s), "\n";
				}
				else
				{
					print $fout join ("\t", $chrom, $chromStart, $chromEnd + 1, $name, $score), "\n";
				}
			}
			$iter++;	
		} #cluster
	}#chrom
}#strand
close ($fout);
close ($fout) if $outBedFile ne '-';

system ("rm -rf $cache");

