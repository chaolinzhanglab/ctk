#!/usr/bin/perl -w

use strict;
use File::Basename;
use Carp;
use Getopt::Long;

use Bed;

my $verbose = 0;

my $prog = basename ($0);
my $rgb = "blue";

my $trackName = "test_track";


GetOptions (
	"name:s"=>\$trackName,
	"col:s"=>\$rgb,
	"v|verbose"=>\$verbose
);

my %col2rgb = (
	red=>	'255,0,0',
	darkred=>'139,0,0',
	green=>	'0,255,0',
	darkgreen=>'0,100,0',
	lightgreen=>'144,238,144',
	blue=>	'0,0,255',
	darkblue=>'0,0,139',
	lightblue=>'173,216,230',
	skyblue=>'135,206,235',
	pink=>'255,192,203',
	purple=>'160,32,240',
	brown=> '165,42,42',
	cyan=>  '0,255,255',
	yellow=> '255,255,0'
);





if (@ARGV != 2)
{
	print "add rgb color to a bed file\n";
	print "$prog [options] <in.bed> <out.bed>\n";
	print " -col [string] : color by name or rgb ([blue]|red|green|...|r,g,b)\n";
	print "                 ", join ("|", sort keys %col2rgb), "\n";

	print " -v            : verbose\n";
	exit (1);
}

my ($inBedFile, $outBedFile) = @ARGV;


$rgb = $col2rgb{$rgb} if exists $col2rgb{$rgb};


my ($fin, $fout);

my $header = "track name= $trackName itemRgb=\"On\"";

open ($fin, "<$inBedFile") || Carp::croak "can not open file $inBedFile to read\n";
open ($fout, ">$outBedFile") || Carp::croak "can not open file $outBedFile to write\n";

my $i = 0;
my $foundHeader = 0;
my $headerPrinted = 0;
while (my $line =<$fin>)
{
	chomp $line;
	next if $line=~/^\s*$/;
	next if $line=~/^\#/;
	if ($line =~/^track name/)
	{
		$header = $line;
		if ($header !~/itemRgb/i)
		{
			$header .= "\titemRgb=\"On\"";
		}
		print $fout $header, "\n";
		$foundHeader = 1;
		next;
	}

	if ($foundHeader && $headerPrinted == 0)
	{
		print $fout $header, "\n";
		$headerPrinted = 1;
	}

	print "$i ...\n" if $verbose && $i % 10000 == 0;
	$i++;

	my $entry = lineToBed ($line);

	if (exists $entry->{"itemRgb"})
	{
		$entry->{"itemRgb"} = $rgb;
	}
	else
	{
		$entry->{"itemRgb"} = $rgb;
		$entry->{"thickStart"} = $entry->{"chromStart"} unless exists $entry->{"thickStart"};
		$entry->{"thickEnd"} = $entry->{"chromEnd"} unless exists $entry->{"thickEnd"};
	}

	print $fout bedToLine ($entry), "\n";
}


close ($fin);
close ($fout);


