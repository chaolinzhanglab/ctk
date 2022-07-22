#!/usr/bin/env perl
#
use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Carp;
use Bed;
use Sequence;

my $prog=basename($0);
my $cmdDir=dirname($0);

my $selectType = "";
my $substType = "";

my $summaryFile = "";
my $verbose = 0;

GetOptions ('t|mutationType:s'=>\$selectType,
			'subt:s'=>\$substType,
			'summary:s'=>\$summaryFile,
			'v'=>\$verbose);


if (@ARGV != 2 && @ARGV != 1)
{
	print "Get specific types of mutations\n";
	print "usage: $prog [options] <tag.mutation.txt> [tag.mutation.type.bed]\n";
	print " <tag.mutation.txt> : input mutation file\n";
	print " [tag.mutation.type.bed] : output mutations of selected type (if not provided, only summary will be printed\n";
	print "options:\n";
	print " -t         [string]: the type of mutations (del|ins|sub)\n";
	print " --subt     [string]: substitution type (e.g., t2c, effective only when using -t sub; all substitutions when not specified)\n";
	print " --summary  [string]: print summary statistics to file\n";
	print " -v                 : verbose\n";
	exit (1);
}


my ($mutationFile, $outBedFile) = @ARGV;
my %summaryHash;

#if ($selectType ne 'del' && $selectType ne 'ins' && $selectType ne 'sub' && $selectType ne 't2c' && $selectType ne 'c2t')
if ($selectType ne 'del' && $selectType ne 'ins' && $selectType ne 'sub' && $selectType ne '')
{
	Carp::croak "Wrong type of mutation \n";
}

my ($substFrom, $substTo);

if ($selectType eq 'sub')
{

	if ($substType ne '')
	{
		my $tmp;
		($substFrom, $tmp, $substTo) =split (//, $substType);
		$substFrom = uc($substFrom);
		$substTo = uc ($substTo);

		$substFrom = 'T' if $substFrom eq 'U';
		$substTo = 'T' if $substTo eq 'U';

		Carp::croak "wrong substitution pattern: $substType\n" unless $substFrom=~/[ACGT]/ && $substTo=~/[ACGT]/ && $substFrom ne $substTo;
	}
}

my ($fin, $fout);

open ($fin, "<$mutationFile") || Carp::croak "cannot open $mutationFile to read\n";

if (defined $outBedFile && $outBedFile ne '')
{
	open ($fout, ">$outBedFile") || Carp::croak "cannot open $outBedFile to write\n";
}

my $iter = 0;
while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;

	print "$iter ...\n" if $verbose && $iter % 100000 == 0;
	$iter++;

	my @cols = split (/\s+/, $line);
	
	my $strand = $cols[5];
	my ($from, $type, $to) = @cols [7..9];


	if ($type eq '-')
	{
		$summaryHash{'del'}++;
		if (defined $outBedFile && $outBedFile ne '')
		{
			print $fout join ("\t", @cols[0..5]), "\n" if $selectType eq 'del';
		}
	}
	elsif ($type eq '+')
	{
		$summaryHash{'ins'}++;
		if (defined $outBedFile && $outBedFile ne '')
        {
			print $fout join ("\t", @cols[0..5]), "\n" if $selectType eq 'ins';
		}
	}
	elsif ($type eq '>')
	{

		if ($strand eq '-')
		{
			$from = uc (revcom ($from));
			$to = uc (revcom ($to));
		}
		else
		{
			$from = uc ($from);
            $to = uc ($to);
		}
		
		$summaryHash{'sub'}{"$from-$to"}++;
		$summaryHash{'sub'}{'total'}++;

		if ($selectType eq 'sub')
		{
			if ($substType ne '')
			{
				next unless $from eq $substFrom && $to eq $substTo;
			}

			if (defined $outBedFile && $outBedFile ne '')
        	{
				print $fout join ("\t", @cols[0..5]), "\n";
			}
		}
	}
}

close ($fin);
if (defined $outBedFile && $outBedFile ne '')
{
	close ($fout);
}

if ($summaryFile ne '')
{
	open ($fout, ">$summaryFile") || Carp::croak "cannot open file $summaryFile to write\n";

	print $fout "Deletion: ", exists $summaryHash{'del'} ? $summaryHash{'del'} : 0, "\n";
	print $fout "Insertion: ", exists $summaryHash{'ins'} ? $summaryHash{'ins'} : 0, "\n";

	if (exists $summaryHash{'sub'})
	{
		print $fout "Substitution: ", $summaryHash{'sub'}{'total'}, "\n";
		foreach my $t (sort keys %{$summaryHash{'sub'}})
		{
			next if $t eq 'total';
			print $fout "$t: ", $summaryHash{'sub'}{$t}, "\n";
		}
	}

	close ($fout);
}


