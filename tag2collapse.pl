#!/usr/bin/perl -w

use strict;
use warnings;

use File::Basename;
use Data::Dumper;
use Getopt::Long;
use Carp;
use POSIX;

use sort 'stable';

use MyConfig;
use Bed 1.02;
use Motif 1.01;
use Common;


my $prog = basename($0);

my $bigFile = 0;	#if yes, we need to use cache
my $minBlockSize = 200000;

#my $keepScore = 0; #keep the original score of representative tag (for collapse mode)


#singleton means no overlap with other regions
my $weight = 0;	#the way to calculate the score in each cluster, without weight means the number of tags, other wise, the sum of scores
my $weightInName = 0; #if --weight is set, weight-in-name tells me to find the weight attached to the sequence name (before degenerate linker if any)

my $randomBarcode = 0; #do not collapse if the barcode sequences are different
my $randomLinker = 0; #the same as above, but obsolete

my $EM = -1; #the threshold for the EM algorithm to report the reads. if 0, do not perform EM, e.g. 50 is a reasonable threshold
my $sequencingErrorModel = 'alignment'; #'em-local'; # or 'alignment', 'fix=0.01', 'em-global'
my $sequencingErrorRate = 0;
my $sequencingErrorFile = "";
my $maxN = 1000; #when exceed, use approximation

my $keepMaxScore = 0; #if not specifiied, the longest one will be kept
my $keepTagName = 0;

my $cache = getDefaultCache ($prog);
my $keepCache = 0;

my $verbose = 0;
my $debug = 0;


my @ARGV0 = @ARGV;

GetOptions (
			'big'=>\$bigFile,
			'weight'=>\$weight,
			'weight-in-name'=>\$weightInName,
			'rl|random-linker'=>\$randomLinker,
			'rb|random-barcode'=>\$randomBarcode,
			'EM:i'=> \$EM,
			'keep-max-score'=>\$keepMaxScore,
			#'keep-score'=>\$keepScore,
			'keep-tag-name'=>\$keepTagName,
			'seq-error-model:s'=>\$sequencingErrorModel,
			'output-seq-error:s'=>\$sequencingErrorFile,
			'c|cache:s'=>\$cache,
			'keep-cache'=>\$keepCache,
			'd|debug'=>\$debug,
			'v|verbose'=>\$verbose
			);


if (@ARGV != 2)
{
	print "collapse tags according to the start position\n";
	
	print "Usage: $prog [options] <in.bed> <out.bed>\n";
	print " <in.bed> : gz file acceptable\n";
	print "[Input options]\n";
	print " -big                        : set when the input file is big\n";
	print " -weight                     : consider the weight of each tag\n";
	print " --weight-in-name            : find weight in name\n";
	print "\n[EM options]\n";
	print " --random-barcode            : random barcode exists, no collapse for different barcodes\n";
	print " -EM        [int]            : EM threshold to infer reliability of each collapsed read (when have random linker, -1=no EM)\n"; 
	print " --seq-error-model  [string] : sequencing error model to use (alignment|[em-local]|em-global|fix=0.01)\n";
	print " --output-seq-error [file]   : output sequencing errors estimated by the EM algorithm\n";
	
	print "\n[Output options]\n";
	print " --keep-max-score            : keep the tag with the most weight (instead of the longest one) as representative\n";
	#print " --keep-score               : keep the original score of representative tags\n";
	print " --keep-tag-name             : do not change tag name (no extra information)\n";

	print "\n[Other options]\n";
	print " -c      [string]            : cache dir ($cache)\n";
	print " --keep-cache                : keep cache when the job is done\n";
	print " -d                          : debug (on|[off])\n";
	print " -v                          : verbose (on|[off])\n";
	exit (1);
}


my ($inBedFile, $outBedFile) = @ARGV;
my $msgio = $outBedFile eq '-' ? *STDERR :  *STDOUT;

$randomBarcode = $randomLinker if ($randomLinker != 0 && $randomBarcode == 0);

print $msgio "CMD = $prog ", join (' ', @ARGV0), "\n"; # if $verbose;


#check consistency of parameters

if ($EM >= 0)
{
	Carp::croak "--random-linker must be set in the EM mode\n" unless $randomBarcode;
	Carp::croak "EM must be <=100\n" if $EM > 100;
}
else
{
	$sequencingErrorModel = '';
}

if ($sequencingErrorModel ne '' || $sequencingErrorFile ne '')
{
	Carp::croak "you must be in the EM mode with --random-linker and --EM {50}\n" unless $randomBarcode && $EM >= 0;
}

if ($sequencingErrorModel=~/^fix/)
{
	if ($sequencingErrorModel=~/^fix=(\S*?)$/)
	{
		$sequencingErrorRate = $1;
		$sequencingErrorModel = 'fix';

		Carp::croak "sequencing error rate must be greater than zero\n" unless $sequencingErrorRate > 0;
	}
	else
	{
		Carp::croak "incorrect format in the --seq-error-model $sequencingErrorModel option\n";
	}
}
elsif ($sequencingErrorModel eq 'alignment')
{
	if ($weight)
	{
		Carp::croak "--weight-in-name must be set and the score column must be mismatches\n" unless $weightInName;
	}
}
elsif ($sequencingErrorModel ne '' && $sequencingErrorModel ne 'em-global' && $sequencingErrorModel ne 'em-local')
{
	Carp::croak "incorrection option --seq-error-model $sequencingErrorModel\n";
}


Carp::croak "$cache already exists\n" if -d $cache;
system ("mkdir $cache");



#read data and split if necessary
my %tagCount;

print "reading tags from $inBedFile ...\n" if $verbose;
if ($bigFile)
{
	my $ret = splitBedFileByChrom ($inBedFile, $cache, v=>$verbose, sort=>1); 
	%tagCount = %$ret;
}
else
{
	#my $tags = readBedFile ($inBedFile, $verbose);

	my $sortedBedFile = $cache . "/". basename ($inBedFile) . "sort";
	sortBedFile ($inBedFile, $sortedBedFile, 0, $cache);
	my $tags = readBedFile ($sortedBedFile, $verbose);
	#we now sort the bed file so that we get the same results with or without big flag (cz 05/05/2018)

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


#estimate sequencing error rate if necessary


my @strands = qw(+ -);

if ($sequencingErrorModel eq 'em-global' || $sequencingErrorModel eq 'alignment')
{
	#estimate the global sequencing error rate
	print "estimating the global sequencing error rate using the $sequencingErrorModel model ...\n";

	if ($sequencingErrorModel eq 'em-global')
	{
		$sequencingErrorRate = 0.02;
		print "setting initial sequencing error rate: $sequencingErrorRate\n" if $verbose;
	}

	my $totalIter = $sequencingErrorModel eq 'alignment' ? 1 : 3;
	for (my $iter = 0; $iter < $totalIter; $iter++)
	{
		my $totalSequencingError = 0;
		my $totalNucleotideCount = 0;
		
		print "\n\niteration $iter ...\n" if $verbose;
		foreach my $s (@strands)
		{
			print "processing strand $s ...\n" if $verbose;
	
			foreach my $chrom (sort keys %tagCount)
			{
				my $tags;
				if ($bigFile)
				{
					my $tmpFile = $tagCount{$chrom}->{'f'};
					print "loading tags on chromsome $chrom from $tmpFile...\n" if $verbose;
					
					my $fin;
					open ($fin, "<$tmpFile") || Carp::croak "cannot open file $tmpFile to read\n";
					
					my $blockId = 0;
					while (my $tags = readNextBedBlock ($fin, 1, 0, minBlockSize=> $minBlockSize))
					{
						next unless @$tags > 0;

						my $ret;
						if ($sequencingErrorModel eq 'em-global')
						{
							#print "update the esitmate of sequencing errors on $chrom ...\n" if $verbose;
							#count the number of sequencing errors in all positions with at least 50 tags
							$ret = updateSeqErrorEM ($tags, $s, $weight, $weightInName, 50, $sequencingErrorRate, $verbose);
						}
						else
						{
							$ret = countSeqError ($tags, $s, $weight, $weightInName, $verbose);
						}

						print $ret->{'E'}, " mismatches in ", $ret->{'NL'}, " nucleotides considered in block $blockId\n" if $verbose;
						
						next unless $ret->{'NL'} > 0;
						$totalSequencingError += $ret->{'E'};
						$totalNucleotideCount += $ret->{'NL'};
						
						$blockId++;
					}

					close ($fin);
				}
				else
				{
					my @tags = @{$tagCount{$chrom}};
					#make a copy rather than use the reference, because the EM algorithm will add stuff
					my $tags = \@tags;
					next unless @$tags > 0;

					my $ret;
					if ($sequencingErrorModel eq 'em-global')
					{
						#print "update the esitmate of sequencing errors on $chrom ...\n" if $verbose;
						#count the number of sequencing errors in all positions with at least 50 tags
						$ret = updateSeqErrorEM ($tags, $s, $weight, $weightInName, 50, $sequencingErrorRate, $verbose);
					}
					else
					{
						$ret = countSeqError ($tags, $s, $weight, $weightInName, $verbose);
					}

					print $ret->{'E'}, " mismatches in ", $ret->{'NL'}, " nucleotides considered on $chrom\n" if $verbose;
			
					next unless $ret->{'NL'} > 0;
					$totalSequencingError += $ret->{'E'};
					$totalNucleotideCount += $ret->{'NL'};
				}
			}#chrom
		}#strand
		print "total error = $totalSequencingError, total nucleotide count=$totalNucleotideCount\n";

		$sequencingErrorRate = $totalSequencingError / $totalNucleotideCount;
		print "updated sequencing error rate: $sequencingErrorRate\n" if $verbose;
	}#iter
}



#now start to collapse tags according to starting coordinates (and random linker if necessary)

print "\n\ncollapsing tags ...\n" if $verbose;

my $fout;
open ($fout, ">$outBedFile") || Carp::croak "can not open file $outBedFile to write\n";

my $fout2;
if ($sequencingErrorFile ne '')
{
	open ($fout2, ">$sequencingErrorFile") || Carp::croak "cannot open file $sequencingErrorFile to write\n";
}


foreach my $s (@strands)
{
	print "processing strand $s ...\n" if $verbose;
	foreach my $chrom (sort keys %tagCount)
	{
		
		#retreive tags
		my $tags;
		if ($bigFile)
		{
			my $tmpFile = $tagCount{$chrom}->{'f'};
			
			my $fin;
			open ($fin, "<$tmpFile") || Carp::croak "cannot open file $tmpFile to read\n";

			my $blockId = 0;
			while (my $tags = readNextBedBlock ($fin, 1, 0, minBlockSize=> $minBlockSize))
			{
				print "processing block # $blockId on $chrom ...\n" if $verbose;
				tag2collapse ($tags, $s, $fout, $fout2);
				$blockId++;
			}
			close ($fin);
		}
		else
		{
			print "processing $chrom ...\n" if $verbose;
			$tags = $tagCount{$chrom};
			tag2collapse ($tags, $s, $fout, $fout2);
		}
	}#chrom
}#strand
close ($fout);

if ($sequencingErrorFile ne '')
{
	close ($fout2);
}

system ("rm -rf $cache") unless $keepCache;





=head2 tag2collapse

output unique tags

=cut

sub tag2collapse
{
	my ($tags, $strand, $fout, $fout2) = @_;
	#tags: tags on a chrom or on a block
	#fout: output tags
	#fout2: output sequencing errors
		
	return unless @$tags > 0;
	
	Carp::croak "No strand specified\n" unless exists $tags->[0]->{"strand"};
	my $chrom = $tags->[0]->{'chrom'};
	my $s = $strand;
	
	if ($weight && $weightInName)
	{
		#put the copy number attached to the name to the score column
		foreach my $t (@$tags)
		{
			my @cols = split (/\#|\-/, $t->{'name'}); #NOTE: fastx_collapser used - to separate read id & the number of copy 
			pop @cols if $randomBarcode;
			$t->{'score'} = pop @cols;
		}
	}

	my $n = @$tags;
	print "$n tags loaded\n" if $verbose;
	
	#clustering
	print "clustering tags by starting coordinates...\n" if $verbose;

	#in each cluster, the reads are sorted by length (descending)
	my $clusters = collapseReads ($tags, $s);

	my $nc = @$clusters;
	print "$nc clusters found\n\n" if $verbose;

	my $iter = 0;
	foreach my $clust (@$clusters)
	{
		print "cluster $iter ...\n" if $debug;
		$iter++;


		#count the number of tags (with copy number considered) in each cluster
		my $score = countTagInCluster ($clust, $weight);			

		if ($weight && $keepMaxScore)
		{
			#the tag with the most number of copies will be ranked first
			my @clustSorted = sort {$b->{'score'} <=> $a->{'score'}} @$clust;
			$clust = \@clustSorted;
		}


		if ($randomBarcode)
		{
			#divide tags according to the random linker
			my %tagsInClust;
			my $linkerIter = 0;
			foreach my $tag (@$clust)
			{
				my $linker = getRandomLinker ($tag->{"name"});
				next if $linker =~/[^ACGT]/i;

				#$linkerIter++ unless exists $tagsInClust{$linker};
				$tagsInClust{$linker}->{'iter'} = $linkerIter++;
				#add linkerIter to allow stable sorting

				push @{$tagsInClust{$linker}->{"tag"}}, $tag;
				my $s = $weight ? $tag->{"score"} : 1;
				$tagsInClust{$linker}->{"score"}+= $s; 
			}
			next unless keys %tagsInClust > 0;

			#EM
			if ($EM >= 0)
			{
				#for the em-local model, the error will be estimated on the fly
				#for all other models, including alignemnt, fix, or em-global, the error rate is already calculated 
				my $EM_summary = {};
				if (keys %tagsInClust <= $maxN)
				{
					$EM_summary = $sequencingErrorModel eq 'em-local' ?
                        mixtureEM (\%tagsInClust) : mixtureEM (\%tagsInClust, $sequencingErrorRate);
				}
				else
				{
					$EM_summary = $sequencingErrorModel eq 'em-local' ?  
						mixtureEMSparse (\%tagsInClust) : mixtureEMSparse (\%tagsInClust, $sequencingErrorRate);
				}
				#output summary information if requested
				print $fout2 join ("\t", $EM_summary->{'k'}, $EM_summary->{'N'}, $EM_summary->{'e'}), "\n" 
				if $sequencingErrorFile ne '';
			

				#output
				my $nuniq = 0; #No of unique reads at this position
				
				my @sortedLinkers = sort {$tagsInClust{$a}{'iter'} <=> $tagsInClust{$b}{'iter'}} keys %tagsInClust;
				
				foreach my $linker (sort {$tagsInClust{$b}->{'error'} <=> $tagsInClust{$a}->{'error'} || $tagsInClust{$b}->{'n'} <=> $tagsInClust{$a}->{'n'}} @sortedLinkers)
				#foreach my $linker (sort {$tagsInClust{$b}->{'error'} <=> $tagsInClust{$a}->{'error'} || $a cmp $b} keys %tagsInClust)
				{	#sort according to confidence

					#print "e=", $tagsInClust{$linker}->{'error'}, "\n";
					my $tagsInClustWithSameLinker = $tagsInClust{$linker};

					#if nothing is significant, output the best one
					#because we believe there is at least one real unique tag
					$tagsInClustWithSameLinker->{"error"} = 100 if $nuniq == 0 && $tagsInClustWithSameLinker->{"error"} < $EM;
						
					#filter if EM was performed
					next if $tagsInClustWithSameLinker->{"error"} < $EM;
					$nuniq++;

					#the representative tag will be the first one, which is the longest one by default, 
					#or the one with most copies with --wieght and --keep-max score were set
					my $representativeTag =  $tagsInClustWithSameLinker->{"tag"}->[0];

					#if (not $keepScore)
					#{
					$representativeTag->{"score"} = sprintf ("%.2f", $tagsInClustWithSameLinker->{"error"});
					#}

					if (not $keepTagName)
					{
						$representativeTag->{"name"} .= "[o=" . $tagsInClustWithSameLinker->{"score"} . "]";
						$representativeTag->{"name"} .= 
							"[n=" .	sprintf ("%.0f", $tagsInClustWithSameLinker->{"n"}). "][e=" . 
							sprintf ("%.0f", $tagsInClustWithSameLinker->{"error"}) . "]" unless $keepTagName;
					}
					print $fout bedToLine ($representativeTag), "\n";
				}
				print "$nuniq unique tags identified in cluster\n" if $debug;
			}
			else
			{
				#No EM
				#output
				foreach my $linker (keys %tagsInClust)
				{
					my $tagsInClustWithSameLinker = $tagsInClust{$linker};
					
					#the representative tag will be the first one, which is the longest one by default, or the one with most copies with --wieght and --keep-max score were set
					my $representativeTag =  $tagsInClustWithSameLinker->{"tag"}->[0];

					$representativeTag->{"score"} = $tagsInClustWithSameLinker->{"score"};
					$representativeTag->{"name"} .= "[o=" . $tagsInClustWithSameLinker->{"score"} . "]" unless $keepTagName;
					print $fout bedToLine ($representativeTag), "\n";
				}
			}
		}
		else
		{	#no barcode
			$clust->[0]->{"score"} = $score; # unless $keepScore;
			print $fout bedToLine ($clust->[0]), "\n";
		}
	} #cluster
}



=head2 mixtureEM

Usage:
my $ret = mixtureEM ($tags); #will try to estimate error rate
or 
my $ret = mixtureEM ($tags, 0.01); #will use the given error rate

return:
	{k=>($#linker+1), #number of unique linkers
	N=>$N, 			  #total number of reads
	L=> $L, 		  #length of linker
	E=> $totalError,  #total number of sequencing errors
	e=> $e};		  #error rate

=cut

sub mixtureEM
{
	my $tags = $_[0];
	
	my $e = 0.02;
	my $fixError = 0;
	if (@_> 1)
	{
		#use the provided error rate and fix it during EM
		$e = $_[1];
		$fixError = 1;
	}

	my @linker = sort keys %$tags;
	my $debugOld = $debug;
	$debug = 1 if @linker > $maxN; #turn on the debug flag for more output when there are a lot of tags and thus it is very slow
	
	Carp::croak "no tags exist for EM algorithm\n" unless @linker > 0;

	my $t = $tags->{$linker[0]}->{'tag'}->[0];
	#Carp::croak Dumper ($t), "\n";
	my $pos = $t->{'strand'} eq '+' ? $t->{'chromStart'} : $t->{'chromEnd'};
	print "chrom=", $t->{'chrom'}, ", strand=", $t->{'strand'}, ", pos=$pos\n" if $debug;

	if (@linker == 1)
	{
		my $linker = $linker[0];
		$tags->{$linker}->{"error"} = "100";
		$tags->{$linker}->{"n"} = $tags->{$linker}->{"score"};
		return {k=>1, N=> $tags->{$linker}->{"score"}, L=> length($linker), e=>0, E=> 0};
	}
	
	my @alpha;		#prior probability of each obs belonging to a component
	my $N = 0;

	print "initializing EM ...\n" if $debug;
	print scalar (@linker), " observations:\n" if $debug;
	for (my $k = 0; $k < @linker; $k++)
	{
		my $linker = $linker[$k];
		$alpha[$k] = $tags->{$linker}->{"score"};
		$N += $alpha[$k];
		print "$k: ", $linker, " [", $tags->{$linker}->{"score"}, "]\n" if $debug;
	}
	
	Carp::croak "N = 0?\n" unless $N > 0;

	for (my $k = 0; $k < @linker; $k++)
	{
		$alpha[$k] /= $N;
	}

	my $L = length ($linker[0]);
	
    my @linkerBase;
    for (my $i = 0; $i < @linker; $i++)
    {
        my $w = $linker[$i];
        $w =~tr/a-z/A-Z/;
        my @b = split (//, $w);
        push @linkerBase, \@b;
    }

	my @alphaPost;	#posterior probability of each obs belonging to a component
	my @dist; #distance (number of mismatches) between linkers
    
	for (my $i = 0; $i < @linker; $i++)
    {
		#preallocate memory of dist and alphaPost
		my @a = ();
		my @b = ();

		$#a = scalar (@linker);
		$#b = scalar (@linker);

		$dist[$i] = \@a;
		$alphaPost[$i] = \@b;
	}

	for (my $i = 0; $i < @linker; $i++)
	{
        $dist[$i][$i] = 0;
        for (my $j = $i+1; $j < @linker; $j++)
        {
            $dist[$i][$j] = countMismatchBase ($linkerBase[$i], $linkerBase[$j]);
            $dist[$j][$i] = $dist[$i][$j];
        }
    }

	print "starting EM iterations ...\n" if $debug;

	my $totalError = 0;
	for (my $iter = 0; $iter < 100; $iter++)
	{
		#E-step
		for (my $i = 0; $i < @linker; $i++)
		{#obs
			my $denominator = 0;
			for (my $k = 0; $k < @linker; $k++)
			{#component
				$alphaPost[$k][$i] = ($e ** $dist[$i][$k]) * ((1 - $e)** ($L-$dist[$i][$k])) * $alpha[$k];
				$denominator += $alphaPost[$k][$i];
			}

			for (my $k = 0; $k < @linker; $k++)
			{
				$alphaPost[$k][$i] /=  $denominator;
			}
		}

		my @alphaOld;

		#update alpha
		for (my $k = 0; $k < @linker; $k++)
		{#compnent
			my $s = 0;
			for (my $i = 0; $i < @linker; $i++)
			{
				$s += $alphaPost[$k][$i] * $tags->{$linker[$i]}->{"score"};
			}
			$alphaOld[$k] = $alpha[$k];
			$alpha[$k] = $s / $N;
		}
	
		#M-step
		#update error
		$totalError = 0;
		for (my $i = 0; $i < @linker; $i++)
		{
			for (my $k = 0; $k < @linker; $k++)
			{
				#$totalError += $tags->{$linker[$i]}->{"score"} * $alphaPost[$k][$i] * ($e ** $dist[$i][$k]) * ((1 - $e)** ($L-$dist[$i][$k])) * $dist[$i][$k];
				$totalError += $dist[$i][$k] * $alphaPost[$k][$i] * $tags->{$linker[$i]}->{"score"};
			}
		}
		
		if (not $fixError)
		{
			$e = $totalError / $N / $L;
			#bound the estimatin of error rate
			$e = max($e, 0.001);
			$e = min($e, 0.02);
		}

		my $improve = 0;
		for (my $k = 0; $k < @linker; $k++)
		{
			$improve += ($alpha[$k] - $alphaOld[$k])**2;
		}

		$improve /= @linker;
		$improve = sqrt ($improve);

		print "iter $iter: estimated sequencing error=$e, improvement=$improve\n" if $debug;
		last if $improve < 1e-6;
	}
	
	print "iteration finished\n" if $debug;

	for (my $k = 0; $k < @linker; $k++)
	{
		my $logP = 0;
		for (my $i = 0; $i < @linker; $i++)
		{
			#Carp::croak "current linker $k=", $linker[$k], ", $i=", $linker[$i], Dumper ($tags), "\n" if $alphaPost[$k][$i] == 1;

			$alphaPost[$k][$i] = 1 - (1e-16)  if $alphaPost[$k][$i] == 1;
			$logP += $tags->{$linker[$i]}->{"score"} * log(1-$alphaPost[$k][$i]) / log(10);
		}
		$tags->{$linker[$k]}->{"error"} = -$logP * 10;
		$tags->{$linker[$k]}->{"n"} = $alpha[$k] * $N;
	}

	for (my $k = 0; $k < @linker; $k++)
	{
		my $linker = $linker[$k];
		print "$k: ", $linker, " [", $tags->{$linker}->{"score"}, "], \trelative abundance=", 
			sprintf ("%.3f", ($alpha[$k] * $N)) ,"\treliability=", 
			sprintf ("%.3f", $tags->{$linker[$k]}->{"error"}) ,"\n" if $debug;
	}
	$debug = $debugOld;
	return {k=>($#linker+1), N=>$N, L=> $L, E=> $totalError, e=> $e};
}


=head2 mixtureEMSparse

#consider only linkers with at most 3 mismatches

Usage:
my $ret = mixtureEMSparse ($tags); #will try to estimate error rate
or 
my $ret = mixtureEMSparse ($tags, 0.01); #will use the given error rate

return:
	{k=>($#linker+1), #number of unique linkers
	N=>$N, 			  #total number of reads
	L=> $L, 		  #length of linker
	E=> $totalError,  #total number of sequencing errors
	e=> $e};		  #error rate

=cut

sub mixtureEMSparse
{
	my $tags = $_[0];
	
	my $e = 0.02;
	my $fixError = 0;
	if (@_> 1)
	{
		#use the provided error rate and fix it during EM
		$e = $_[1];
		$fixError = 1;
	}

	my @linker = sort keys %$tags;
	my $debugOld = $debug;
	$debug = 1 if @linker > $maxN; #turn on the debug flag for more output when there are a lot of tags and thus it is very slow
	
	Carp::croak "no tags exist for EM algorithm\n" unless @linker > 0;

	my $t = $tags->{$linker[0]}->{'tag'}->[0];
	#Carp::croak Dumper ($t), "\n";
	my $pos = $t->{'strand'} eq '+' ? $t->{'chromStart'} : $t->{'chromEnd'};
	print "chrom=", $t->{'chrom'}, ", strand=", $t->{'strand'}, ", pos=$pos\n" if $debug;

	if (@linker == 1)
	{
		my $linker = $linker[0];
		$tags->{$linker}->{"error"} = "100" + rand (1e-5);
		$tags->{$linker}->{"n"} = $tags->{$linker}->{"score"};
		return {k=>1, N=> $tags->{$linker}->{"score"}, L=> length($linker), e=>0, E=> 0};
	}
	
	my @alpha;		#prior probability of each obs belonging to a component
	my $N = 0;

	print "initializing EM ...\n" if $debug;
	print scalar (@linker), " observations:\n" if $debug;
	
	my $maxCopy = 0; #the barcode with the max copy number
	for (my $k = 0; $k < @linker; $k++)
	{
		my $linker = $linker[$k];
		$alpha[$k] = $tags->{$linker}->{"score"};
		$N += $alpha[$k];
		$maxCopy = $alpha[$k] if $maxCopy < $alpha[$k];

		print "$k: ", $linker, " [", $tags->{$linker}->{"score"}, "]\n" if $debug;
	}
	
	Carp::croak "N = 0?\n" unless $N > 0;

	for (my $k = 0; $k < @linker; $k++)
	{
		$alpha[$k] /= $N;
	}

	my $L = length ($linker[0]);
	
    my @linkerBase;
    for (my $i = 0; $i < @linker; $i++)
    {
        my $w = $linker[$i];
        $w =~tr/a-z/A-Z/;
        my @b = split (//, $w);
        push @linkerBase, \@b;
    }

	my @alphaPost;	#posterior probability of each obs belonging to a component
	my @dist; #distance (number of mismatches) between linkers
    
	my $maxDist = 3;

	#see if we can use smaller maxDist to save time, change in 04/23/2016
	if (@linker > 2000)
	{
		#too low and too much memory with larger maxDist 
		$maxDist = 1;
	}
	elsif (@linker > $maxN)
	{
		$maxDist = ceil(log($maxCopy)/log(1/$e));
		$maxDist = 3 if $maxDist > 3;
		#$maxDist = 2 if $maxDist < 2;
	}
	print "e=$e, maxCopy=$maxCopy, maxDist = $maxDist\n" if $debug;

	for (my $i = 0; $i < @linker; $i++)
	{
        $dist[$i]{$i} = 0;
        for (my $j = $i+1; $j < @linker; $j++)
        {
			my $mc = countMismatchBase ($linkerBase[$i], $linkerBase[$j], $maxDist+1);
			#we count at most $maxDist+1 differences to save time
			if ($mc <= $maxDist)
			{
				#we assume at most $maxDist mismaches can occur
            	$dist[$i]{$j} = $mc;
            	$dist[$j]{$i} = $dist[$i]{$j};
			}
        }
    }

	print "starting EM iterations ...\n" if $debug;

	my $totalError = 0;
	for (my $iter = 0; $iter < 100; $iter++)
	{
		#E-step
		for (my $i = 0; $i < @linker; $i++)
		{#obs
			my $denominator = 0;
			foreach my $k (keys %{$dist[$i]})
			{#component
				$alphaPost[$k]{$i} = ($e ** $dist[$i]{$k}) * ((1 - $e)** ($L-$dist[$i]{$k})) * $alpha[$k];
				$denominator += $alphaPost[$k]{$i};
			}

			foreach my $k (keys %{$dist[$i]})
			{	
				#0 otherwise
				$alphaPost[$k]{$i} /=  $denominator;
			}
		}

		my @alphaOld;

		#update alpha
		for (my $k = 0; $k < @linker; $k++)
		{#compnent
			my $s = 0;
			
			foreach my $i (keys %{$dist[$k]})
			{
				$s += $alphaPost[$k]{$i} * $tags->{$linker[$i]}->{"score"};
			}
			$alphaOld[$k] = $alpha[$k];
			$alpha[$k] = $s / $N;
		}
	
		#M-step
		#update error
		$totalError = 0;
		for (my $i = 0; $i < @linker; $i++)
		{
			foreach my $k (keys %{$dist[$i]})
			{
				#$totalError += $tags->{$linker[$i]}->{"score"} * $alphaPost[$k][$i] * ($e ** $dist[$i][$k]) * ((1 - $e)** ($L-$dist[$i][$k])) * $dist[$i][$k];
				$totalError += $dist[$i]{$k} * $alphaPost[$k]{$i} * $tags->{$linker[$i]}->{"score"};
			}
		}
		
		if (not $fixError)
		{
			$e = $totalError / $N / $L;
			#bound the estimatin of error rate
			$e = max($e, 0.001);
			$e = min($e, 0.02);
		}

		my $improve = 0;
		for (my $k = 0; $k < @linker; $k++)
		{
			$improve += ($alpha[$k] - $alphaOld[$k])**2;
		}

		$improve /= @linker;
		$improve = sqrt ($improve);

		print "iter $iter: estimated sequencing error=$e, improvement=$improve\n" if $debug;
		last if $improve < 1e-6;
	}
	
	print "iteration finished\n" if $debug;

	for (my $k = 0; $k < @linker; $k++)
	{
		my $logP = 0;
		
		foreach my $i (keys %{$dist[$k]})
		{
			#Carp::croak "current linker $k=", $linker[$k], ", $i=", $linker[$i], Dumper ($tags), "\n" if $alphaPost[$k][$i] == 1;
			$alphaPost[$k]{$i} = 1 - (1e-16)  if $alphaPost[$k]{$i} == 1;
			$logP += $tags->{$linker[$i]}->{"score"} * log(1-$alphaPost[$k]{$i}) / log(10);
		}
		$tags->{$linker[$k]}->{"error"} = -$logP * 10;
		$tags->{$linker[$k]}->{"n"} = $alpha[$k] * $N;
	}

	for (my $k = 0; $k < @linker; $k++)
	{
		my $linker = $linker[$k];
		print "$k: ", $linker, " [", $tags->{$linker}->{"score"}, "], \trelative abundance=", 
			sprintf ("%.3f", ($alpha[$k] * $N)) ,"\treliability=", 
			sprintf ("%.3f", $tags->{$linker[$k]}->{"error"}) ,"\n" if $debug;
	}
	$debug = $debugOld;
	return {k=>($#linker+1), N=>$N, L=> $L, E=> $totalError, e=> $e};
}







sub countTagInCluster
{
	my ($clust, $weight) = @_;

	my $score = $weight ? 0 : @$clust;
	map {$score += $_->{"score"}} @$clust if $weight;
	return $score;
}


#update the estimate of sequencing errors for all tags on a chromosome
#
#assume we have degenerate linker in the name and the number of mismatch is in the score
#and the copy number is attached to the name
sub countSeqError
{
	my ($tags, $strand, $weight, $weightInName, $verbose) = @_;

	if ($weight)
	{
		Carp::croak "weight must be attached in seq name\n" unless $weightInName;
	}

	my $E = 0;
	my $NL = 0;
	foreach my $t (@$tags)
	{
		next if $strand ne 'b' && $t->{'strand'} ne $strand;

		my $mismatch = $t->{'score'};

		#a data consistency check, better than nothing
		Carp::croak "incorrect number of mismatches, tag=", Dumper ($t), "\n" 
		unless $mismatch >= 0 && $mismatch < $t->{'chromEnd'} - $t->{'chromStart'} + 1;
		
		my $copyNum = 1;

		if ($weight && $weightInName)
		{	
			my @cols = split (/\#|\-/, $t->{'name'}); #NOTE: fastx_collapser used - to separate read id & the number of copy 
			Carp::croak "no copy number or degenerate linker? ", Dumper ($t), "\n" if @cols < 3;
			pop @cols;
			$copyNum = pop @cols;
		}
		$E += $mismatch * $copyNum;

		#update on Feb 3 2016 to accommdate exon junction reads
		my $s = $t->{'chromEnd'} - $t->{'chromStart'} + 1;
		$s = sum($t->{'blockSizes'}) if exists $t->{'blockSizes'};
		$NL += $s * $copyNum;
	}
	return {E=>$E, NL=>$NL};
}


sub updateSeqErrorEM
{
	#NOTE: tags are assumed to be on the same chromosome and only those on the given strand are considered
	my ($tags, $strand, $weight, $weightInName, $minTagInCluster, $currErrorRate, $verbose) = @_;
	if ($weight && $weightInName)
	{
		foreach my $t (@$tags)
		{
			my @cols = split (/\#|\-/, $t->{'name'}); #NOTE: fastx_collapser used - to separate read id & the number of copy 
			pop @cols; #assume we have random linker here
			$t->{'score'} = pop @cols;
		}
	}

	my $n = @$tags;
	
	Carp::croak "No strand specified\n" unless exists $tags->[0]->{"strand"};

	print "clustering tags ...\n" if $debug;

	#cluster reads according to genomic starting coordinates
	#in each cluster, the reads are sorted by length (descending)
	my $clusters = collapseReads ($tags, $strand);

	my $nc = @$clusters;
	print "\n\n$nc clusters found\n" if $debug;

	my $iter = 0;
	my %ret;

	my $N = 0; #total number of tags
	my $E = 0; #total number of mutations
	my $L = 0; #linker lenegth
	foreach my $clust (@$clusters)
	{
		my $score = countTagInCluster ($clust, $weight);
		next unless $score >= $minTagInCluster;

		#further group tags according to the degenerate linker
		#and record the score of each group
		my %tagsInClust;
		foreach my $tag (@$clust)
		{
			my $linker = getRandomLinker ($tag->{"name"});
			next if $linker =~/[^ACGT]/i;

			$L = length ($linker);
			push @{$tagsInClust{$linker}->{"tag"}}, $tag;

			my $s = $weight ? $tag->{"score"} : 1;
			$tagsInClust{$linker}->{"score"}+= $s; 
		}
				
		my $EM_summary = {};
		if (keys %tagsInClust <= $maxN)
		{
			$EM_summary = mixtureEM (\%tagsInClust, $currErrorRate);
		}
		else
		{
			$EM_summary = mixtureEMSparse (\%tagsInClust, $currErrorRate);
		}

		if ($EM_summary->{'N'} >= $minTagInCluster)
		{
			$E += $EM_summary->{'E'};
			$N += $EM_summary->{'N'};
		}
	}
	return {E=> $E, N=> $N, NL=>$N * $L};
}

sub getRandomLinker
{
	my $tagName = $_[0];
	my @cols = split (/\#|\-/, $tagName);
	Carp::croak "no random linker is found in $tagName\n" if @cols < 2;
	my $linker = pop @cols;
	return $linker;
}

