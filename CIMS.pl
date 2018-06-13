#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use File::Basename;

use MyConfig;
use Bed;


my $prog = basename ($0);

my $verbose = 0;

my $mutationSize = 1;
my $numIter = 5;
my $trackPosition = 0;
my $mutationPosSummaryFile = ""; #output position of mutations
my $big = 0;

my $noSparseCorrect = 0;
my $FDR = 1;
my $mkr = 0.00;

my $cache = getDefaultCache ($prog);
my $keepCache = 0;

my @ARGV0 = @ARGV;

GetOptions ("n:i"=>\$numIter,
		'w:i'=>\$mutationSize,
		"c:s"=>\$cache,
		'keep-cache'=>\$keepCache,
		'p'=>\$trackPosition,
		'outp:s'=>\$mutationPosSummaryFile,
		'big'=>\$big,
		'no-sparse-correct'=>\$noSparseCorrect,
		'FDR:f'=>\$FDR,
		'mkr:f'=>\$mkr,
		"v"=>\$verbose);



if (@ARGV != 3)
{
	print "identify significant CIMS\n";
	print "Usage: $prog [options] <tag.bed> <mutation.bed> <out.txt>\n";
	print "Note: the 5th column of <mutation.bed> must have the position of the mutation relative to the chromStart of the tag\n";
	print " -big               : big file\n";
	print " -w     [int]       : mutation size ($mutationSize)\n";
	print " -n     [int]       : number of iterations for permutation ($numIter)\n";
	print " -p                 : track mutation position relative to read start\n";
	print " --outp [file]      : output mutation position summary\n";
	print " --no-sparse-correct: no sparcity correction\n";
	print " -FDR   [float]     : threshold of FDR ($FDR)\n";
	print " -mkr   [float]     : threshold of m-over-k-ratio ($mkr)\n";
	print " -c     [dir]       : cache dir ($cache)\n";
	print " --keep-cache       : keep cache when the job is done\n";
	print " -v                 : verbose\n";
	exit (1);
}


my ($tagBedFile, $mutationBedFile, $outFile) = @ARGV;

my $cmdDir = dirname ($0);


print "CMD = $prog ", join (" ", @ARGV0), "\n" if $verbose;

my $verboseFlag = $verbose ? '-v' : '';
my $bigFlag = $big ? '-big' :'';

system ("mkdir $cache"); 


print "clean mutations ...\n" if $verbose;
my $cmd = "awk \'{if(\$3-\$2==$mutationSize) {print \$0}}\' $mutationBedFile > $cache/mutation.clean.bed";
system ($cmd);

$mutationBedFile = "$cache/mutation.clean.bed";

$cmd = "wc -l $mutationBedFile | awk '{print \$1}'";
my $n = `$cmd`;
$n=~/^(\d+)/;
$n = $1;
Carp::croak "no mutaton of size $mutationSize found\n" unless $n > 0;

print "$n mutations of size $mutationSize found\n";

my $mutationPositionFile = "$cache/mutation.pos.txt";
if ($trackPosition)
{
	#the original position is always relative to chromStart (0 based), even for reads on the negative strand
	#the awk convert the position relative to read start (on the sense strand)
	
	print "generating CIMS position file $mutationPositionFile...\n" if $verbose;
	
	my $tmpMutationBedFile = "$cache/mutation.bed";
	my $cmd = "sort -k 4,4 $mutationBedFile > $tmpMutationBedFile";
	print "$cmd\n";	
	system ($cmd);
	$mutationBedFile = $tmpMutationBedFile;	 

	#TODO:
	#need to perform data integrity check to make sure the mutation coordinates are within the tags


	my $mutationIdFile = "$cache/mutation.id";

	$cmd = "cut -f 4 $mutationBedFile | awk -F \"[\" '{print \$1}' > $mutationIdFile";
	#my $cmd = "cut -f 4 $mutationBedFile | awk -F \"[\" '{print \$1}' | sort > $mutationIdFile";
	
	print $cmd, "\n" if $verbose;
	system ($cmd);

	my $tmpFile = "$cache/tmp";
	$cmd = "python $cmdDir/joinWrapper.py $tagBedFile $mutationIdFile 4 1 N $tmpFile";
	#$cmd = "perl ~/scripts/selectRow.pl -q 3 $tagBedFile $mutationIdFile > $tmpFile";
	print "$cmd\n" if $verbose;
	my $ret = system ($cmd);
	print "CMD=$cmd failed: $?\n" if $ret != 0;

	my $mutationNum = `wc -l $mutationIdFile | awk '{print \$1}'`; chomp $mutationNum;
	my $tagNum = `wc -l $tmpFile | awk '{print \$1}'`; chomp $tmpFile;

	if ($mutationNum != $tagNum)
	{
		Carp::croak "Inconsistency detected between the tag file and the mutation file. Make sure the NAME column of the two files use the same IDs\n";
	}
	$cmd = "paste $mutationBedFile $tmpFile | awk '{if(\$6 ==\"+\") {print \$5;} else {print \$9-\$8\-\$5-$mutationSize}}' > $mutationPositionFile";
	#my $cmd = "awk '{if(\$6 ==\"+\") {print \$5;} else {print \$3-\$2-\$5-1}}' $mutationBedFile > $mutationPositionFile";
	system ($cmd);

	if ($mutationPosSummaryFile ne '')
	{
		$cmd = "cut -f 2 $mutationPositionFile | sort -n | uniq -c | awk '{print \$2\"\\t\"\$1}'> $mutationPosSummaryFile";
		$ret = system ($cmd);
		print "CMD=$cmd failed: $?\n" if $ret != 0;
	}
}

print "clustering mutation sites ...\n" if $verbose;
my $mutationClusterFile = "$cache/mutation.cluster.bed";
$cmd = "perl $cmdDir/tag2cluster.pl  -s  -maxgap \"-$mutationSize\" $bigFlag $verboseFlag $mutationBedFile $mutationClusterFile";
$cmd .= " -big" if $big;
system ($cmd);

#
print "counting the total tag number for each mutation site ...\n" if $verbose;
my $mutationTagCountFile = "$cache/mutation.tagcount.bed";
$cmd = "perl $cmdDir/tag2profile.pl -c $cache/tmp_mut_count -region $mutationClusterFile -ss -of bed $bigFlag $verboseFlag $tagBedFile $mutationTagCountFile";
my $ret = system ($cmd);
print "CMD=$cmd failed: $?\n" if $ret != 0;

#
print "reading frequency of mutations and tag number for each mutation site\n" if $verbose;
my $mutationClusters = readBedFile ($mutationClusterFile, $verbose);  # m
my $mutationTagCount = readBedFile ($mutationTagCountFile, $verbose); # k

my $totalMutations = 0;
my %mutationClusterOverallHash;
my %mutationClusterHash;

for (my $i = 0; $i < @$mutationClusters; $i++)
{
	my $c = $mutationClusters->[$i];
	$c->{'freq'} = $c->{'score'};
	$totalMutations += $c->{'freq'};

	$c->{'tag'} = $mutationTagCount->[$i]->{'score'};

	$mutationClusterHash{$c->{'tag'}}{$c->{'freq'}}->{'count'}++;
	$mutationClusterOverallHash{$c->{'freq'}}->{'count'}++;
}


print "find cummulative distributions ...\n" if $verbose;

my $s = 0;
foreach my $freq (sort {$b <=> $a} keys %mutationClusterOverallHash)
{
	$mutationClusterOverallHash{$freq}->{'cumm'} = $s + $mutationClusterOverallHash{$freq}->{'count'};
	$s = $mutationClusterOverallHash{$freq}->{'cumm'};
}

foreach my $tagNum (sort {$b <=> $a} keys %mutationClusterHash)
{
	my $mutationClusterHashWithTagNum = $mutationClusterHash{$tagNum};
	my $s = 0;
	foreach my $freq (sort {$b<=>$a} keys %$mutationClusterHashWithTagNum)
	{
		$mutationClusterHashWithTagNum->{$freq}->{'cumm'} = $s + $mutationClusterHashWithTagNum->{$freq}->{'count'};
		$s = $mutationClusterHashWithTagNum->{$freq}->{'cumm'};
	}
}



#my $numIter = 1;

srand (0);

print "start simulations, $totalMutations mutations in total  ...\n" if $verbose;

my %randomMutationClusterHash;
my %randomMutationClusterOverallHash;
for (my $i = 0; $i < $numIter; $i++)
{
	print "Iteration $i ...\n" if $verbose;

	#
	print "generating random CIMS sites ...\n" if $verbose;
	my $randomMutationBedFile = "$cache/mutation.random.$i.bed";
	my $trackPosFlag = $trackPosition ? " -p $mutationPositionFile " : "";
	my $rand_seed = 0;

	#generate a non-zero seed for each iteration
	$rand_seed = rand () + 1;

	#$cmd = "perl ~/scripts/simulateCLIPmismatch.pl -n $totalMutations $trackPosFlag $verboseFlag -srand $rand_seed $tagBedFile $randomMutationBedFile";
	$cmd = "perl $cmdDir/simulateCIMS.pl -n $totalMutations -w $mutationSize $trackPosFlag -srand $rand_seed $tagBedFile $randomMutationBedFile";
	
	print $cmd, "\n" if $verbose;
	my $ret = system ($cmd);

	Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;

	#
	print "clustering mutation sites in $randomMutationBedFile ...\n" if $verbose;
	my $randomMutationClusterFile = "$cache/mutation.random.$i.cluster.bed";
	#$cmd = "perl ~/scripts/tag2cluster.pl  -s  -maxgap \"-1\" $verboseFlag $randomMutationBedFile $randomMutationClusterFile";
	$cmd = "perl $cmdDir/tag2cluster.pl  -s  -maxgap \"-$mutationSize\" $randomMutationBedFile $randomMutationClusterFile";
	
	$cmd .= " -big" if $big;
	$ret = system ($cmd);

	Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;

	#
	print "counting the total tag number for each clustered mutation site in $randomMutationClusterFile ...\n" if $verbose;
	my $randomMutationTagCountFile = "$cache/mutation.random.$i.tagcount.bed";
	#$cmd = "perl ~/scripts/tag2profile.pl -region $randomMutationClusterFile -ss -of bed $verboseFlag $tagBedFile $randomMutationTagCountFile";
	$cmd = "perl $cmdDir/tag2profile.pl -c $cache/tmp_mut_count_$i -region $randomMutationClusterFile -ss -of bed $tagBedFile $randomMutationTagCountFile";
	
	$cmd .= " -big" if $big;
	$ret = system ($cmd);
	Carp::croak "CMD=$cmd failed: $?\n" if $ret != 0;

	print "reading frequency of mutations and tag number for each mutation site\n" if $verbose;

	#my $randomCIMSTagCount = AnnotationIO::readBedFile ($randomCIMSTagCountFile, $verbose);
	#my $randomCIMSClusters = AnnotationIO::readBedFile ($randomCIMSClusterFile, $verbose);

	my $fin1;
	my $fin2;
	open ($fin1, "<$randomMutationTagCountFile") || Carp::croak "cannot open file $randomMutationTagCountFile to read\n";
	open ($fin2, "<$randomMutationClusterFile") || Carp::croak "cannot open file $randomMutationClusterFile to read\n";

	my $i = 0;
	while (my $line = <$fin1>)
	{
		#my $line = <$fin1>;
		chomp $line;
		my @cols = split ("\t", $line);
		my $k = $cols[4];

		$line = <$fin2>;
		chomp $line;
		@cols = split ("\t", $line);
		my $m = $cols[4];

		$randomMutationClusterHash{$k}{$m}->{'count'}++;
		$randomMutationClusterOverallHash{$m}->{'count'}++;

		$i++;
	}
	close ($fin1);
	close ($fin2);
}

print "find cummulative distributions ...\n" if $verbose;
$s = 0;
foreach my $freq (sort {$b <=> $a} keys %randomMutationClusterOverallHash)
{
	$randomMutationClusterOverallHash{$freq}->{'cumm'} = $s + $randomMutationClusterOverallHash{$freq}->{'count'};
	$s = $randomMutationClusterOverallHash{$freq}->{'cumm'};
}

foreach my $tagNum (sort {$b <=> $a} keys %randomMutationClusterHash)
{
	my $randomMutationClusterHashWithTagNum = $randomMutationClusterHash{$tagNum};
	my $s = 0;
	foreach my $freq (sort {$b<=>$a} keys %$randomMutationClusterHashWithTagNum)
	{
		$randomMutationClusterHashWithTagNum->{$freq}->{'cumm'} = $s + $randomMutationClusterHashWithTagNum->{$freq}->{'count'};
		$s = $randomMutationClusterHashWithTagNum->{$freq}->{'cumm'};
	}
}

my (%qHash, %qHash2);
foreach my $c (@$mutationClusters)
{
	my $tagNum = $c->{'tag'};
	my $freq = $c->{'freq'};

	my $count = exists $mutationClusterHash{$tagNum} && exists $mutationClusterHash{$tagNum}{$freq} ? $mutationClusterHash{$tagNum}{$freq}->{'cumm'} : 0;
	my $countOverall = exists $mutationClusterOverallHash{$freq} ? $mutationClusterOverallHash{$freq}->{'cumm'} : 0;

	my $countRandom = exists $randomMutationClusterHash{$tagNum} && exists $randomMutationClusterHash{$tagNum}{$freq} ? $randomMutationClusterHash{$tagNum}{$freq}->{'cumm'} : 0;
	my $countOverallRandom = exists $randomMutationClusterOverallHash{$freq} ? $randomMutationClusterOverallHash{$freq}->{'cumm'} : 0;

	$countRandom /= $numIter;
	$countOverallRandom /= $numIter;

	my $q = $countRandom / $count;
	$q = 1 if $q >1;
	my $qOverall = $countOverallRandom / $countOverall;
	$qOverall = 1 if $qOverall > 1;

	$c->{'name'} .= "[k=$tagNum][m=$freq]";
	$c->{'q'} = $q;

	$qHash{$freq}{$tagNum} = \$q;
	$qHash2{$tagNum}{$freq} = \$q;
	#print $fout join ("\t", bedToLine ($c), $tagNum, $freq, $q, $count), "\n" #, $qOverall, $countOverall), "\n";
}

if ($noSparseCorrect == 0)
{
	print "correcting FDRs for sparse sampling ...\n" if $verbose;

	for (my $iter = 0; $iter < 1; $iter++)
	{
		foreach my $m (sort {$a <=> $b} keys %qHash)
		{
			my @ks = sort {$a <=> $b} keys %{$qHash{$m}};
			for (my $i = 1; $i < @ks; $i++)
			{
				my $k1 = $ks[$i-1];
				my $k2 = $ks[$i];
				${$qHash{$m}{$k2}} = ${$qHash{$m}{$k1}} if ${$qHash{$m}{$k2}} < ${$qHash{$m}{$k1}};
			}
		}
=extra
		foreach my $k (sort {$a <=> $b} keys %qHash2)
		{
			my @ms = sort {$b <=> $a} keys %{$qHash2{$k}};
			for (my $i = 1; $i < @ms; $i++)
			{
				my $m1 = $ms[$i-1];
				my $m2 = $ms[$i];
				${$qHash2{$k}{$m2}} = ${$qHash2{$k}{$m1}} if ${$qHash2{$k}{$m2}} < ${$qHash2{$k}{$m1}};
			}
		}
=cut
	}
}



print "writing output to $outFile ...\n" if $verbose;
my $fout;

open ($fout, ">$outFile") || Carp::croak "can not open file $outFile to write\n";
print $fout "#", join ("\t", "chrom", "chromStart", "chromEnd", "name", "score", "strand", "tagNumber(k)", "mutationFreq(m)", "FDR", "count(>=m,k)"), "\n"; #"overallFDR", "overallCount(>=m)"), "\n";

foreach my $c (@$mutationClusters)
{
	my $tagNum = $c->{'tag'};
	my $freq = $c->{'freq'};
	my $q = ${$qHash{$freq}{$tagNum}};
	my $count = exists $mutationClusterHash{$tagNum} && exists $mutationClusterHash{$tagNum}{$freq} ? $mutationClusterHash{$tagNum}{$freq}->{'cumm'} : 0;
	print $fout join ("\t", bedToLine ($c), $tagNum, $freq, $q, $count), "\n" if $q <= $FDR && $freq / $tagNum >= $mkr; #, $qOverall, $countOverall), "\n";
}

close ($fout);

system ("rm -rf $cache") unless $keepCache;
#select(STDOUT);




