#!/usr/bin/perl -w
#
use strict;
use Getopt::Long;
use File::Basename;
use Carp;


my $prog = basename ($0);

my $colId = 0;
my $filterColId = 0;
my $ignoreCase = 0;
my $reverse = 0;
my $verbose = 0;


GetOptions ('q|query-column-id:i'=>\$colId,
		'f:i'=>\$filterColId,
		'i|ignore-case'=>\$ignoreCase,
		'r|reverse'=>\$reverse,
		'v|verbose'=>\$verbose
);


if (@ARGV != 2)
{
	print "remove selected rows from a file\n";
	print "Usage: $prog [options] <input> <filter>\n";
	print " <input>  : use - for STDIN\n";
	print "OPTIONS:\n";
	print " -q [int] : query column id (zero-based) (default=$colId)\n";
	print " -f [int] : filter column id (zero-based) (default=$filterColId)\n";
	print " -i       : ignore case (default=off)\n";
	print " -r       : reverse mode\n";
	print " -v       : verbose\n";
	exit (0);
}

my ($inputFile, $filterFile) = @ARGV;
print STDERR "reading $filterFile ...\n" if $verbose;

my $fin;
open ($fin, "<$filterFile") || die "can not open file $filterFile to read\n";

my %filter;
my $line;
my $i = 0;
while ($line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;

#	if ($line =~/^(.*?)\t/)
#	{
#		$line = $1;
#	}

	my @cols = split (/\t/, $line);

	Carp::croak "line=$line does not have enough columns\n" unless @cols > $filterColId;

	my $line = $cols[$filterColId];

	$i++;
	print STDERR "$i ...\n" if $verbose && $i % 100000 == 0;

	$line=~tr/a-z/A-Z/ if $ignoreCase;
	
	$filter{$line} = 1;
}
close ($fin);


print STDERR "reading $inputFile ...\n" if $verbose;

if ($inputFile eq '-')
{
	$fin = *STDIN;
}
else
{
	open ($fin, "<$inputFile") || die "can not open file $inputFile to read\n";
}

$i = 0;
while ($line = <$fin>)
{
	chomp $line;
	my @cols = split (/\s+/, $line);
	next if $#cols < $colId;
	
	$i++;
	
	print STDERR "$i ...\n" if $verbose && $i % 100000 == 0;
	my $key = $cols[$colId];
	$key =~tr/a-z/A-Z/ if $ignoreCase;
	if ($reverse)
	{
		print "$line\n" if exists $filter{$key};
	}
	else
	{
		print "$line\n" unless exists $filter{$key};
	}
}
close ($fin) if $inputFile ne '-';




