#!/usr/bin/perl -w
#
use strict;
use warnings;

use Carp;
use Getopt::Long;
use File::Basename;

my $prog = basename ($0);

my $printNoMatch = 0;
my $printNoMatchText = "No match found";
my $printSingleLine = 0;
my $printOneMatch = 0;
my $colId = 0;
my $filterColId = 0;
my $ignoreCase = 0;
my $delimitor = "||";
my $withHeader = 0;

GetOptions ('q|query-column-id:i'=>\$colId,
		'f:i'=>\$filterColId,
		'i|ignore-case'=>\$ignoreCase,
		'p|print-no-match'=>\$printNoMatch,
		'pt|print-no-match-text:s'=>\$printNoMatchText,
		's|print-single-line'=>\$printSingleLine,
		'd|delimitor:s'=>\$delimitor,
		'h|with-header'=>\$withHeader,
		'ss'=>\$printOneMatch
);

if (@ARGV != 2)
{
	print "Select rows from a file\n";
	print "Usage: $prog [options] <input-file> <filter-file>\n";
	print " <input-file>: use - for stdin\n";
	print "OPTIONS:\n";
	print " -h       : the file has header to be included\n";
	print " -q [int] : query column id (zero-based) (default=$colId)\n";
	print " -f [int] : filter column id (zero-based) (default=$filterColId)\n";
	print " -i       : ignore case (default=off)\n";
	print " -p       : print the key without matches (default=off)\n";
	print " -pt      : the text to print when no match exist ($printNoMatchText)\n";
	print " -s       : print in single line if there are multiple matches (default=off)\n";
	print " -d       : delimitor ($delimitor) when -s is specified\n";
	print " -ss      : print only one match if there are multiple matches (default=off)\n"; 
	exit (1);
}

my ($input, $filter) = @ARGV;

my $fin;

if ($input eq '-')
{
    $fin = *STDIN;
}
else
{
	open ($fin, "<$input") || Carp::croak "can not open file $input to read\n";
}

my %input;

my $header = "";

if ($withHeader)
{
	while ($header = <$fin>)
	{
		chomp $header;
		last unless $header =~/^\s*$/;
	}
}

my $line;
while ($line = <$fin>)
{
	chomp $line;
	next if $line=~/^\s*$/;
	my @cols = split (/\t/, $line);
	#print "numcol = ", $#cols, "\n";
	#print join ("\n", @cols), "\n\n";
	next if ($colId >$#cols);
	#Carp::croak "column id =$colId, do not have so many columns.\n" if ($colId > $#cols);
	my $key = $cols[$colId];
	$key=~tr/a-z/A-Z/ if ($ignoreCase);
	#print $key, "\n";
	push @{$input{$key}}, $line;	#in case multiple lines have the same key
}
close ($fin);

print $header, "\n" if $withHeader;

open ($fin, "<$filter") || Carp::croak "can not open file $filter to read\n";
while ($line = <$fin>)
{
	chomp $line;
	if ($line=~/^\s*$/)
	{
		if ($printNoMatch)
		{
			print "\n";
			next;
		}
		else
		{
			next;
		}
	}
	
	my @cols = split (/\t/, $line);

	if (@cols <= $filterColId || $cols[$filterColId] eq '')
	{
		if ($printNoMatch)
		{
			print "\n";
			next;
		}
		else
		{
			next;
		}
	
	}
	my $key = $cols[$filterColId];
	#print $key, "\n";
	$key=~tr/a-z/A-Z/ if ($ignoreCase);
	if (exists $input{$key})
	{
		if ($printOneMatch)
		{
			my $outLine = $input{$key}->[0];
			print $outLine, "\n";
		}
		elsif ($printSingleLine)
		{
			my $outLine = joinLines ($input{$key}, $colId, $delimitor); #join ($delimitor, @{$input{$key}});
			print $outLine, "\n";
		}
		else
		{
			foreach (@{$input{$key}})
			{
				print $_, "\n";
			}
		}
	}
	elsif ($printNoMatch)
	{
		my $prefix="\t"x $colId;
		print $prefix . "$key\t$printNoMatchText\n";
	}
}
close ($fin) if $input ne '-';

sub joinLines 
{
	my ($lines, $colId, $delimitor) = @_;
	my @outItems;

	foreach my $line (@$lines)
	{
		my @cols = split (/\t/, $line);
		for (my $i = 0; $i < @cols; $i++)
		{
			if ($i < @outItems && $i != $colId)
			{
				if ($outItems[$i] eq '')
				{
					$outItems[$i] = $cols[$i];
				}
				elsif ($cols[$i] ne '')
				{
					$outItems[$i] = join ($delimitor, $outItems[$i], $cols[$i]);
				}
				#elsif (cols[$i] eq '') do nothing
			}
			else
			{
				$outItems[$i] = $cols[$i];
			}
		}	
	}
	return join ("\t", @outItems);
}

