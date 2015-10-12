#!/usr/bin/perl -w

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Carp;
use Data::Dumper;

use Bed;


my ($leftExt, $rightExt) = (0, 0);
my $filter = ""; #5utr, 3utr, coding
my $tolPartialCoding = 0;
my $outExonFile = "";
my $outIntronFile = "";
my $internal = 0;
my $terminal = 0;
my $terminal3 = 0;
my $terminal5 = 0;
my $completeCodonOnly = 0;
my $noExonIntronId = 0;
my $verbose = 0;
my $itype = 'all';
my $ilenThreshold = 240;
my $keepLncRNA = 0;
my $keepSingleNt = 0;

my $allowDupName = 0;

my $prog = basename ($0);

GetOptions ('l|left:i'=>\$leftExt,
		'r|right:i'=>\$rightExt,
		'f|filter:s'=>\$filter,
		'partial'=>\$tolPartialCoding,
		'itype:s' =>\$itype,
		'ilen:i'=>\$ilenThreshold,
		'oi:s'=>\$outIntronFile,
		'oe:s'=>\$outExonFile,
		'internal'=>\$internal,
		'terminal'=>\$terminal,
		'terminal3' =>\$terminal3,
		'terminal5' =>\$terminal5,
		'cc'=>\$completeCodonOnly,
		'nid'=>\$noExonIntronId,
		'keep-lncrna'=>\$keepLncRNA,
		'keep-1-nt'=>\$keepSingleNt,
		'allow-duplicate-name'=>\$allowDupName,
		'v|verbose'=>\$verbose);

if (@ARGV != 1)
{
	print "get exons or introns from a gene bed file\n";
	print "Usage: $prog [options] <gene.bed>\n";
	print "OPTIONS:\n";
	print " -f     [string]        : filter of region (5utr|coding|3utr)\n";
	print " -itype [string]        : (flanking) intron type ([all]|short|long|)\n";
	print " -ilen     [int]        : threshold of (flanking) intron length ($ilenThreshold)\n";
	print " -partial               : tolerate partial coding when filtering (off)\n";
	print " -internal              : only internal exons (off)\n";
	print " -terminal              : only terminal exons (off)\n";
	print " -terminal3             : only the 3' terminal exon (off)\n";
	print " -terminal5             : only the 5' terminal exon (off)\n";
	print " -l        [int]        : extention on the left ($leftExt)\n";
	print " -r        [int]        : extention on the right ($rightExt)\n";
	print " -cc                    : complete codon only (effective only when filtering with coding) (off)\n";
	print " -oi    [string]        : output intron file\n";
	print " -oe    [string]        : output exon file\n";
	print " --allow-duplicate-name : allow duplicate gene name\n";
	print " -nid                   : no exon or intron id should be appended to transcript accession\n";
	print " --keep-lncrna          : treat lncRNA as utrs\n";
	print " --keep-1-nt            : keep 1-nt exons\n";
	print " -v                     : verbose\n"; 
	exit (0);
}

Carp::croak "internal and terminal can not be specified at the same time\n" if $internal + $terminal + $terminal3 + $terminal5 > 1;

Carp::croak "output nothing?\n" if ($outExonFile eq '' && $outIntronFile eq '');

Carp::croak "$outExonFile already exists\n" if ($outExonFile ne '' && -f $outExonFile);
Carp::croak "$outIntronFile already exists\n" if ($outIntronFile ne '' && -f $outIntronFile);

warn "extension will not be applied to exons\n" if $filter eq 'coding' && ($leftExt != 0 && $rightExt != 0) && $verbose;


my $inBedFile = $ARGV[0];



#by default, we will check duplicate gene names, and if there are any, we will append suffix to differentiate them
print "examining duplicate names ...\n" if $verbose && $allowDupName == 0;

my $cmd = "awk '{print \$4}' $inBedFile| sort | uniq -c | awk '{if(\$1>1) {print \$2}}'";

my @namesWithDuplicates = $allowDupName ? () : `$cmd`;
chomp @namesWithDuplicates;

#Carp::croak join ("\n", @namesWithDuplicates), "\n";

#now we store all ids with duplicates
my %geneHash = map {$_=>0} @namesWithDuplicates;

#Carp::croak Dumper (\%geneHash), "\n";

my $n = keys %geneHash;
print "$n genes with duplicates\n" if $verbose && $allowDupName == 0;


print "extracting exons/introns from gene file $inBedFile ...\n" if $verbose;


my $fin;

open ($fin, "<$inBedFile");

my ($foutExon, $foutIntron);

if ($outExonFile ne '')
{
	open ($foutExon, ">$outExonFile") || Carp::croak "can not open file $outExonFile to write\n";
}

if ($outIntronFile ne '')
{
	open ($foutIntron, ">$outIntronFile") || Carp::croak "can not open file $outIntronFile to write\n";
}

my $iter = 0;
while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\s*$/;

	next if $line =~/^track\sname/;

	print "$iter ...\n" if $verbose && $iter % 100000 == 0;
	$iter++;

	my $g = lineToBed ($line);

	my $name = $g->{"name"};

	if (exists $geneHash{$name})
	{
		#print "found duplicates for $name\n";
		#now we add the suffix if we know it has duplicates in the file
		$geneHash{$name}++;
		$g->{"name"} = $name . "@" . ($geneHash{$name}-1);
	}

	$name = $g->{"name"};

	determineCompleteCodon ($g) if $filter eq 'coding' && $completeCodonOnly;

	my $chrom = $g->{"chrom"};
	my $chromStart = $g->{"chromStart"};
	my $chromEnd = $g->{"chromEnd"};
	#my $name = $g->{"name"};

	Carp::croak "no strand information\n" unless exists $g->{"strand"};
	
	$g = bedToFull ($g) unless exists $g->{'blockStarts'} && $g->{"blockSizes"};

	#Carp::croak "no block information: $line\n" unless exists $g->{"blockStarts"};
	#Carp::croak "no block information\n" unless exists $g->{"blockSizes"};
	
	my $score = exists $g->{"score"} ? $g->{"score"} : 0;

	my $strand = $g->{"strand"};
	my $blockSizes = $g->{"blockSizes"};
	my $blockStarts = $g->{"blockStarts"};
	my $blockNum = @$blockSizes;
	
	#determine the filter
	my $filterStart = $g->{"chromStart"};
	my $filterEnd = $g->{"chromEnd"};
	
	if (not $keepLncRNA)
	{
		next if ($filter && $g->{'thickEnd'} - $g->{'thickStart'} <= 0);
		#noncoding transcript
	}

	if ($filter eq '5utr')
	{
		if ($strand eq '+')
		{
			$filterEnd = $g->{"thickStart"} - 1;
		}
		else
		{
			$filterStart = $g->{"thickEnd"} + 1;
		}
	}
	elsif ($filter eq '3utr')
	{
		if ($strand eq '+')
		{
			$filterStart = $g->{"thickEnd"} + 1;
		}
		else
		{
			$filterEnd = $g->{"thickStart"} - 1;
		}
	}
	elsif ($filter eq 'coding')
	{
		$filterStart = $g->{"thickStart"};
		$filterEnd = $g->{"thickEnd"};
	}
	elsif ($filter ne '')
	{
		Carp::croak "incorrect filter: $filter\n";
	}
	
	#output exons
	if ($outExonFile ne '')
	{
		
		for (my $j = 0; $j < $blockNum; $j++)
		{
			my $i = $j;
			if ($strand eq '-')
			{
				$i = $blockNum -1 - $j;
			}
		
			my $exonStart = $chromStart + $blockStarts->[$i];
			my $exonEnd = $chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1;
			
			if ($internal)
			{
				next unless $i > 0 && $i < $blockNum - 1;
			}

			if ($terminal)
			{
				next unless ($i == 0 || $i == $blockNum - 1);
			}
			elsif ($terminal3)
			{
				next unless $j == $blockNum - 1;
			}
			elsif ($terminal5)
			{
				next unless $j == 0;
			}

			if ($tolPartialCoding)#tolerate partial coding exon
			{
				next if ($exonEnd < $filterStart || $exonStart > $filterEnd);
			}
			else
			{
				next unless ($exonStart >= $filterStart && $exonEnd <= $filterEnd);
			}

			if ($completeCodonOnly && $filter eq 'coding')
			{
				next unless (exists $g->{"chopStart"} && exists $g->{"chopEnd"});

				next unless $g->{"chopStart"}->[$i] >= 0 && $g->{"chopEnd"}->[$i] >= 0;
				$exonStart += $g->{"chopStart"}->[$i];
				$exonEnd -= $g->{"chopEnd"}->[$i];
			}
			
			$exonStart = $filterStart if $exonStart < $filterStart;
			$exonEnd = $filterEnd if $exonEnd > $filterEnd;

			if ($completeCodonOnly == 0 || $filter ne 'coding')
			{
				#we do not want to apply extension if we want to get complete codons
				if ($strand eq '+')
				{
					$exonStart -= $leftExt;
					$exonEnd += $rightExt;
				}
				else
				{
					$exonStart -= $rightExt;
					$exonEnd += $leftExt;
				}
			}
			
			if (not $keepSingleNt)
			{
				next if ($exonStart >= $exonEnd); #to be consistent, we don't want exons with only 1nt
			}

			#filtering according to flanking intron length
			my ($uilen, $dilen) = (0, 0); #upstream and downstream intron length
			if ($ilenThreshold > 0 && $itype ne 'all')
			{
				$uilen = $blockStarts->[$i] - ($blockStarts->[$i-1] + $blockSizes->[$i-1] - 1) if ($i > 0);
				$dilen = $blockStarts->[$i+1] - ($blockStarts->[$i] + $blockSizes->[$i] - 1) if ($i < $blockNum -1);
				($uilen, $dilen) = ($dilen, $uilen) if ($strand eq '-');
				
					#this does not matter, just to make more acurate

				if ($itype eq 'long')
				{
					next unless $uilen >= $ilenThreshold && $dilen >= $ilenThreshold;
				}
				elsif ($itype eq 'short')
				{
					next unless $uilen <= $ilenThreshold && $dilen <= $ilenThreshold;
				}
				else
				{
					Carp::croak "illegal value for itype: $itype\n";
				}
			}

			if ($noExonIntronId)
			{
				print $foutExon join ("\t", $chrom, $exonStart, $exonEnd +1, $name, $score, $strand), "\n";
			}
			else
			{
				print $foutExon join ("\t", $chrom, $exonStart, $exonEnd +1, $name ."_$j", $score, $strand), "\n"; #exons and introns are numbered from 5' to 3'
			}
		}
	}
	
	#output introns
	if ($outIntronFile ne '')
	{
		for (my $j = 0; $j < $blockNum; $j++)
		{
			my $i = $j;
			if ($strand eq '-')
			{
				$i = $blockNum -2 - $j;
			}
		
			next if $i+1 > $blockNum - 1;
	
			my $intronStart = $chromStart + $blockStarts->[$i] + $blockSizes->[$i];
			my $intronEnd = $chromStart + $blockStarts->[$i+1] - 1;
			
			next unless ($intronStart >= $filterStart && $intronEnd <= $filterEnd);
			
			my $intronLen = $intronEnd - $intronStart + 1;
			if ($strand eq '+')
			{
				$intronStart -= $leftExt;
				$intronEnd += $rightExt;
			}
			else
			{
				$intronStart -= $rightExt;
				$intronEnd += $leftExt;
			}

			next if ($intronStart > $intronEnd);

			#filtering according to intron length
			if ($itype eq 'long')
			{
				next unless $intronLen >= $ilenThreshold;
			}
			elsif ($itype eq 'short')
			{
				next unless $intronLen <= $ilenThreshold;
			}
			elsif ($itype ne 'all')
			{
				Carp::croak "illegal value for itype: $itype\n";
			}

			if ($noExonIntronId)
			{
				print $foutIntron join ("\t", $chrom, $intronStart, $intronEnd +1, $name, $score, $strand), "\n";
			}
			else
			{
				print $foutIntron join ("\t", $chrom, $intronStart, $intronEnd +1, $name ."_$j", $score, $strand), "\n"; #exon and introns are numbered from 5' to 3'
			}
		}
	}
}


if ($outExonFile ne '')
{
	close ($foutExon);
}

if ($outIntronFile ne '')
{
	close ($foutIntron);
}

#this subroutine has some bugs (when the CDS is in a single exon of a multi-exon gene)
#need to be replaced by the one in Bed.pm (determineExonReadingFrame)

sub determineCompleteCodon
{

	my $g = $_[0];

	my $chrom = $g->{"chrom"};
	my $chromStart = $g->{"chromStart"};
	my $chromEnd = $g->{"chromEnd"};
	my $name = $g->{"name"};
	
	#next if $name ne 'NM_207396';
	
	Carp::croak "no strand information\n" unless exists $g->{"strand"};
	Carp::croak "no block information\n" unless exists $g->{"blockStarts"};
	Carp::croak "no block information\n" unless exists $g->{"blockSizes"};
	
	my $strand = $g->{"strand"};
	my $blockSizes = $g->{"blockSizes"};
	my $blockStarts = $g->{"blockStarts"};
	my $blockNum = @$blockSizes;
	
	my $thickStart = $g->{"thickStart"};
	my $thickEnd = $g->{"thickEnd"};

	for (my $i = 0; $i < $blockNum; $i++)
	{
		$g->{"chopStart"}->[$i] = -1;
		$g->{"chopEnd"}->[$i] = -1;
	}

	if ($thickStart == $chromStart || $thickEnd == $chromEnd)
	{
		warn "$name: ORF is probably incomplete\n" if $verbose;
		next;
	}
	
	#intronless gene
	if ($blockNum == 1)
	{
		my $s = $g->{"blockSizes"}->[0] - ($thickStart - $chromStart) - ($chromEnd - $thickEnd);
		if ($s - int (($s+0.5)/3) * 3 != 0)
		{
			warn "$name: ORF does not have a length divided by 3\n" if $verbose;
			next;
		}
		
		$g->{"chopStart"}->[0] = $thickStart - $chromStart;
		$g->{"chopEnd"}->[0] = $chromEnd - $thickEnd;
		next;
	}

	#print Dumper ($g), "\n";
	#determine the length of cds
	my $cdsLen = 0;

	for (my $i = 0; $i < $blockNum; $i++)
	{
		if($chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 < $thickStart)
		{
			#5'/3' UTR exon
			#next;
		}
		elsif ($chromStart + $blockStarts->[$i] < $thickStart 
				&& $chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 >= $thickStart)
		{
			#the first (partial) coding exon
			$cdsLen = $chromStart + $blockStarts->[$i] + $blockSizes->[$i] - $thickStart;
		}
		elsif ($chromStart + $blockStarts->[$i] >= $thickStart && 
				$chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 <= $thickEnd)
		{
			#complete coding exon
			$cdsLen += $blockSizes->[$i];
		}
		elsif ($chromStart + $blockStarts->[$i] >= $thickStart &&
				$chromStart + $blockStarts->[$i] <= $thickEnd &&
				$chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 > $thickEnd)
		{
			#the last partial coding exon
			$cdsLen += $thickEnd - ($chromStart + $blockStarts->[$i]) + 1;
			last;
		}
		elsif ($chromStart + $blockStarts->[$i] > $thickEnd)
		{
			#3'/5' UTR exon
		}
		else
		{
			Carp::croak "$name: block=$i, something wrong not handled properly...\n";
		}
	}

	#print "cds Len = $cdsLen\n";
	if ($cdsLen % 3 != 0)
	{
		warn "$name: ORF does not have a length divided by 3\n" if $verbose;
		next;
	}

	#determine the complete codons
	my $currCDSLen = 0;
	
	for (my $i = 0; $i < $blockNum; $i++)
	{
		#5'/3'utr exon
		if ($chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 < $thickStart)
		{
			$g->{"chopStart"}->[$i] = -1;
			$g->{"chopEnd"}->[$i] = -1;
		}
		elsif ($chromStart + $blockStarts->[$i] < $thickStart 
				&& $chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 >= $thickStart)
		{
			#the first (partial) coding exon
			$g->{"chopStart"}->[$i] = $thickStart - ($chromStart + $blockStarts->[$i]);
			$currCDSLen = $chromStart + $blockStarts->[$i] + $blockSizes->[$i] - $thickStart;
			$g->{"chopEnd"}->[$i] = $currCDSLen - int(($currCDSLen+0.5)/3) * 3;
			
			my $s = $g->{"blockSizes"}->[$i] - $g->{"chopStart"}->[$i] - $g->{"chopEnd"}->[$i];

			my $exonStart = $chromStart + $blockStarts->[$i]; # + $g->{"chopStart"}->[$i];
			my $exonEnd = $chromStart + $blockStarts->[$i] + $blockSizes->[$i]; # - $g->{"chopEnd"}->[$i];
		
			Carp::croak "$chrom:$exonStart-$exonEnd: $name, block = $i, CDS len = $s, can not be divided by 3\n" 
			unless $s - int (($s+0.5)/3)*3 == 0;
		}
		elsif ($chromStart + $blockStarts->[$i] >= $thickStart && 
				$chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 <= $thickEnd)
		{
			#complete coding exon
			if ($g->{"chopEnd"}->[$i-1] >= 0)
			{
				$g->{"chopStart"}->[$i] = Common::mod (3 - $g->{"chopEnd"}->[$i-1], 3);
			}
			else
			{
				$g->{"chopStart"}->[$i] = 0;
			}
			
			$currCDSLen += $blockSizes->[$i];
			$g->{"chopEnd"}->[$i] = $currCDSLen - int(($currCDSLen+0.5)/3) * 3;
			my $s = $g->{"blockSizes"}->[$i] - $g->{"chopStart"}->[$i] - $g->{"chopEnd"}->[$i];
	
			my $exonStart = $chromStart + $blockStarts->[$i]; # + $g->{"chopStart"}->[$i];
			my $exonEnd = $chromStart + $blockStarts->[$i] + $blockSizes->[$i]; # - $g->{"chopEnd"}->[$i];
			Carp::croak "$chrom:$exonStart-$exonEnd: $name, block = $i, CDS len = $s, can not be divided by 3\n" 
			unless $s - int (($s+0.5)/3)*3 == 0;
	
		}
		elsif ($chromStart + $blockStarts->[$i] >= $thickStart &&
				$chromStart + $blockStarts->[$i] <= $thickEnd &&
				$chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 > $thickEnd)
		{
			#the last partial coding exon
			if ($g->{"chopEnd"}->[$i-1] >= 0)
			{
				$g->{"chopStart"}->[$i] = Common::mod (3 - $g->{"chopEnd"}->[$i-1], 3);
			}
			else
			{
				$g->{"chopStart"}->[$i] = 0;
			}
			
			$g->{"chopEnd"}->[$i] = $chromStart + $blockStarts->[$i] + $blockSizes->[$i] - 1 - $thickEnd;

			my $s = $g->{"blockSizes"}->[$i] - $g->{"chopStart"}->[$i] - $g->{"chopEnd"}->[$i];
			
			my $exonStart = $chromStart + $blockStarts->[$i]; # + $g->{"chopStart"}->[$i];
			my $exonEnd = $chromStart + $blockStarts->[$i] + $blockSizes->[$i]; # - $g->{"chopEnd"}->[$i];
		
			Carp::croak "$chrom:$exonStart-$exonEnd: $name, block = $i, CDS len = $s, can not be divided by 3\n" 
			unless $s - int (($s+0.5)/3)*3 == 0;
		}
		elsif ($chromStart + $blockStarts->[$i] > $thickEnd)
		{
			#3'/5' UTR exon
			$g->{"chopStart"}->[$i] = -1;
			$g->{"chopEnd"}->[$i] = -1;
		}
		else
		{
			Carp::croak "$name: block=$i, something wrong not handled properly...\n";
		}
	}
}


