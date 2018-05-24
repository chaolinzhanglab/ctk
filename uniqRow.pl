#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use File::Basename;

use Common;

my $prog = basename ($0);

my $criteria = "random"; #max_num, min_num, max_text, min_text
my $groupColIdx = 0;	#the column used to group rows
my $compareColIdx = 1;	#the column used to compare rows
my $noExtraCol = 0;
my $oneBasedIndex = 0; #col ids are 1-based, for compatibility with galaxy
my $verbose = 0;

my $withHeader = 0;
my $naString = "NA";

my @ARGV0 = @ARGV;

GetOptions ('c:s'=>\$criteria,
	'h'=>\$withHeader,
	'id:i'=> \$groupColIdx,
	'value:i'=>\$compareColIdx,
	'one-based-index'=>\$oneBasedIndex,
	'no-extra-col'=>\$noExtraCol,
	'na-str:s'=>\$naString,
	'v'=> \$verbose);


if (@ARGV != 2)
{
	print "select uniq rows\n";
	print "Usage: $prog [options] <in.txt> <out.txt>\n";
	print " <in.txt>         : use \"-\" for stdin\n";
	print " <out.txt>        : use \"-\" for stdout\n";
	print "OPTIONS:\n";
	print " -h               : with header line to be included\n";
	print " -c     [string]  : criteria to sort ([random]|max_num|max_num_abs|min_num|min_num_abs|max_text|min_text|sum|mean|rowsum|rowmean|count)\n";
	print " -id    [int]     : id column (from 0) used to group rows ($groupColIdx)\n";
	print " -value [int]     : value column (from 0) used to compare rows ($compareColIdx)\n";
	print " --one-based-index: 1-based column index\n";
	print " --no-extra-col   : do not print extra columns other than the id and value columns\n";
	print " --na-str [string]: NA string ($naString)\n";
	print " -v               : verbose\n";
	exit (1);
}


my ($inFile, $outFile) = @ARGV;


my $msgio = $outFile ne '-' ? *STDOUT : *STDERR;

print $msgio "CMD= $prog ", join (" ", @ARGV0), "\n";
if ($oneBasedIndex)
{
	$groupColIdx--;
	$compareColIdx--;
}

my $fin;

my %rowHash;

if ($inFile eq '-')
{
	$fin = *STDIN;
}
else
{
	open ($fin, "<$inFile") || Carp::croak "can not open file $inFile to read\n";
}


my $header = "";

if ($withHeader)
{
	$header = <$fin>;
	chomp $header;
}

my $nrow = 0;

print $msgio "reading data from $inFile ...\n" if $verbose;
my $i = 0;
while (my $line = <$fin>)
{
	chomp $line;
	next if $line=~/^\s*$/;

	print $msgio "$i ...\n" if $i % 100000 == 0 && $verbose;

	$i++;
	my @cols = split (/\t/, $line);
	my $n = @cols - 1;

	Carp::croak "illegal group col id: $groupColIdx > $n\n" if $groupColIdx > $n;
	Carp::croak "illegal value col id: $compareColIdx > $n\n" if $compareColIdx > $n;

	$nrow++;
	my $groupId = $cols [$groupColIdx];
	my $value = $cols [$compareColIdx];

	push @{$rowHash{$groupId}}, {value=> $value, row=>$line};
}
close ($fin) if $inFile ne '-';

my $ngroup = keys %rowHash;

print $msgio "$nrow rows, $ngroup uniq rows loaded\n" if $verbose;


my $fout;

if ($outFile eq '-')
{
	$fout = *STDOUT;
}
else
{
	open ($fout, ">$outFile") || Carp::croak "can not open file $outFile to write\n";
}

print $msgio "dumping unique rows ...\n" if $verbose;


print $fout $header, "\n" if $withHeader;

$i = 0;
foreach my $groupId (sort keys %rowHash)
{
	my $rows = $rowHash{$groupId};
	my $nr = @$rows;

	my $row = "";

	print $msgio "$i ...\n" if $i % 10000 == 0 && $verbose;
	$i++;

	if ($criteria eq 'random')
	{
		my $idx = randSeq (0, $nr);
		$idx = $idx->[0];
		$row = $rows->[$idx];
	}
	elsif ($criteria eq 'max_num')
	{
		my @rows = sort {$b->{"value"} <=> $a->{"value"}} @$rows;
		$row = $rows[0];
	}
	elsif ($criteria eq 'max_num_abs')
	{
		my @rows = sort {ABS($b->{"value"}) <=> ABS($a->{"value"})} @$rows;
		$row = $rows[0];
	}
	elsif ($criteria eq 'min_num')
	{
		my @rows = sort {$a->{"value"} <=> $b->{"value"}} @$rows;
		$row = $rows[0];
	}
	elsif ($criteria eq 'min_num_abs')
	{
		my @rows = sort {ABS($a->{"value"}) <=> ABS($b->{"value"})} @$rows;
		$row = $rows[0];
	}
	elsif ($criteria eq 'max_text')
	{
		my @rows = sort {$b->{"value"} cmp $a->{"value"}} @$rows;
		$row = $rows[0];
	}
	elsif ($criteria eq 'min_text')
	{
		my @rows = sort {$a->{"value"} cmp $b->{"value"}} @$rows;
		$row = $rows[0];
	}
	elsif ($criteria eq 'sum')
	{
		my $score = 0;
		map {$score += $_->{"value"}} @$rows;

		#for (my $i = 0; $i < @$rows; $i++)
		#{
		#	$score += $rows->[$i]->{"value"};
		#}
		$row = $rows->[0];
		$row->{"value"} = $score;

		my $line = $row->{"row"};
		my @cols = split (/\t/, $line);
		$cols[$compareColIdx] = $score;
		$line = join ("\t", @cols);
		$row->{"row"} = $line;
		#$row = $rows->[0];
	}
	elsif ($criteria eq 'count')
	{
		my $score = @$rows;
		$row = $rows->[0];
		$row->{"value"} = $score;

		my $line = $row->{"row"};
		my @cols = split (/\t/, $line);
		$cols[$compareColIdx] = $score;
		$line = join ("\t", @cols);
		$row->{"row"} = $line;
		
	}
	elsif ($criteria eq 'mean')
	{
		my $score = 0;
		map {$score += $_->{"value"}} @$rows;
		#for (my $i = 0; $i < @$rows; $i++)
		#{
		#	$score += $rows->[$i]->{"value"};
		#}

		$score /= @$rows if @$rows > 0;

		$row = $rows->[0];
		$row->{"value"} = $score;

		my $line = $row->{"row"};
		my @cols = split (/\t/, $line);
		$cols[$compareColIdx] = $score;
		$line = join ("\t", @cols);
		$row->{"row"} = $line;

	}
	elsif ($criteria eq 'rowsum' || $criteria eq 'rowmean')
	{
		$row = $rows->[0];

		my (@cols, @counts); # = split ("\t", $row->{'row'});

		for (my $i = 0; $i < @$rows; $i++)
		{
			my @cols2 = split ("\t", $rows->[$i]->{'row'});
			$cols[$groupColIdx] = $cols2[$groupColIdx];

			for (my $j = 0; $j < @cols2; $j++)
			{
				$counts[$j] += 0;
				next if $j == $groupColIdx || $cols2[$j] eq $naString;

				$cols[$j]+= $cols2[$j];
				$counts[$j] += 1;
			}
		}
		
		my $n = @$rows;
		if ($criteria eq 'rowmean')
		{
			for (my $j = 0; $j < @counts; $j++)
            {
                next if $j == $groupColIdx;
				
                $cols[$j] = $counts[$j] > 0 ? $cols[$j]/$counts[$j] : 'NA';
            }
		}
		
		$row->{"row"} = join ("\t", @cols);
	}
	
	if ($noExtraCol)
	{
		print $fout join ("\t", $groupId, $row->{'value'}), "\n";
	}
	else
	{
		print $fout $row->{"row"}, "\n";
	}
}
close ($fout) if $outFile ne '-';
