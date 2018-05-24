#!/usr/bin/perl -w

use strict;
use warnings;

use Data::Dumper;
use Carp;
use File::Basename;
use Getopt::Long;

use Bed;
use Sam;

my $prog = basename ($0);
my $mapQual = 0;
my $minLen = 0;
my $splitDel = 0;
my $mutationFile = "";
my $verbose = 0;
my $minAnchor = 1; #mutation is considered only when a match 
my $minAnchorIndel = 5;

#my $useRNAStrand = 0; # use the strand of the RNA instead of the read
my $inDelInScore = 0;

GetOptions (
	"map-qual:i"=>\$mapQual,
	"min-len:i"=>\$minLen,
	"split-del"=>\$splitDel,
	"indel-to-end:i"=>\$minAnchorIndel,
	"indel-in-score"=>\$inDelInScore,
	"mutation-file:s"=>\$mutationFile,
	"v|verbose"=>\$verbose);

if (@ARGV != 2)
{
	print STDERR "parse alignment in SAM format to BED format (both tags and mutations)\n";
	print STDERR "Usage: $prog [options] <in.sam> <out.bed>\n";
	print STDERR " <in.sam>  : gzip compressed input file with .gz extension is allowed\n";
	print STDERR " <out.bed> : output bed file of tag alignment\n";
	print STDERR " You can also use - to specify STDIN for input or STDOUT for output\n\n";
	print STDERR "options:\n";
	print STDERR " --mutation-file [string]: file name to save mutations\n";
	print STDERR " --map-qual      [int]   : MAPQ score (e.g. to keep only uniquely mapped reads)\n";
	print STDERR " --min-len       [int]   : minimal read length to report\n";
	print STDERR " --indel-to-end  [int]   : nt from indel to end ($minAnchorIndel)\n";
	print STDERR " --split-del             : split oligo deletion into single nucleotides\n";
	print STDERR " --indel-in-score        : count indels in as mismatches reported in the score column\n";
	print STDERR " -v                      : verbose\n";
	exit (1);
}

my ($inSAMFile, $outBedFile) = @ARGV;
my $msgio = $outBedFile eq '-' ? *STDERR :  *STDOUT;
my ($fin, $fout, $fout2);

if ( $inSAMFile eq "-")
{
    $fin = *STDIN;
}
else
{
	if ($inSAMFile =~/\.gz$/)
	{
		open ($fin, "gunzip -c $inSAMFile | ") || Carp::croak "cannot open file $inSAMFile to read\n";
	}
	else
	{
    	open ($fin, "<$inSAMFile") || Carp::croak "cannot open file $inSAMFile to read\n";
	}
}
if ( $outBedFile eq "-")
{
     $fout = *STDOUT;
}
else
{
    open ($fout, ">$outBedFile") || Carp::croak "cannot open file $outBedFile to write\n";
}

if ($mutationFile ne '')
{
	open ($fout2, ">$mutationFile") || Carp::croak "cannot open file $mutationFile to write\n";
}


my $i = 0;

while (my $line = <$fin>)
{
	chomp $line;

	next if $line=~/^\s*$/;
	next if $line=~/^\@/;

	print $msgio "$i ...\n" if $verbose && $i % 50000 == 0;
	$i++;

	my $sam = lineToSam ($line); # a very lightweight parser
	next unless length($sam->{'SEQ'}) >= $minLen;
	next unless $sam->{"MAPQ"} >=  $mapQual;

	my $flagInfo = decodeSamFlag ($sam->{"FLAG"});
	my $ret = samToBedAndMutation ($sam);
	next unless $ret;

	my $bed = $ret->{'tag'};
	my $mutations = $ret->{'mutations'};

	print $fout bedToLine ($bed), "\n" unless $flagInfo->{'query_nomap'};
	if ($mutationFile ne '')
	{
		foreach my $m (@$mutations)
		{
			my $type = '>'; #substitution
			if ($m->{'type'} eq 'I')
			{
				$type = "+";
			}
			elsif ($m->{'type'} eq 'D')
			{
				$type = '-';
			}
			#output compatible with the previous parser
			if ($type eq '-' && $splitDel)
			{
				my $size = $m->{'chromEnd'} - $m->{'chromStart'} + 1;
				if ($m->{'flag'} >= 2)
				{
					for (my $i = 0; $i < $size; $i++)
					{
						print $fout2 join("\t", $m->{'chrom'}, $m->{'chromStart'} + $i, $m->{'chromStart'} + $i + 1, $m->{'name'}, $m->{'score'} + $i, $m->{'strand'},
							$m->{'score'} + $i, uc(substr ($m->{'refBase'}, $i, 1)), $type, '.', $bed->{'matchStart'}), "\n";
					}
				}
			}
			else
			{
				print $fout2 bedToLine ($m), "\t", join ("\t", $m->{'score'}, uc($m->{'refBase'}), $type, uc($m->{'altBase'}), $bed->{"matchStart"}), "\n" if $m->{'flag'} >= 2;
			}
		}
	}	
}

close ($fin) if $inSAMFile ne '-';
close ($fout) if $outBedFile ne '-';
close ($fout2) if $mutationFile ne '';

sub parseMD 
{
	my ($tagStr, $CIGAR, $SEQ) = @_;

	#insertions are not reflected in MD and they must be extracted from CIGAR
	#so we get CIGAR althought it is not used for now

	my %MDitems = (D=>{}, X=>[]); #substitutions and deletions, but no insertions in MD tag
	if ($tagStr=~/MD\:Z\:(\S+)/g)
	{
		my $MDStr = $1;
		my @token = split (/(\d+)/, $MDStr);
		
		#print join ("\t", "t=", @token), "\n";
	
		next if @token < 2; #no mismatches
		
		#track if there is a good match on the left of each mutation
		my $goodMatchLeft = 0;
		my $goodMatchLeftIndel = 0;
		for (my $i = 0; $i < @token; $i++)
		{
			my $t = $token[$i];
			$token[$i] = {t=>$t};

			if ($t =~/^\d+$/)
			{
				$goodMatchLeft = 1 if $t >= $minAnchor;
				$goodMatchLeftIndel = 1 if $t >= $minAnchorIndel;
			}
			$token[$i]->{'f'} = ($t=~/^\^/) ? $goodMatchLeftIndel : $goodMatchLeft;
		}

		#track if there is a good match on the right of each mutation
		my $goodMatchRight = 0;
		my $goodMatchRightIndel = 0;
		for (my $j = 0; $j < @token; $j++)
		{
			my $i = @token - $j - 1;
			my $t = $token[$i];
			if ($t->{'t'} =~/^\d+$/)
			{
				$goodMatchRight = 1 if $t->{'t'} >= $minAnchor;
				$goodMatchRightIndel = 1 if $t->{'t'} >= $minAnchorIndel;
			}
			$token[$i]{'f'} += ($t->{'t'}=~/^\^/) ? $goodMatchRightIndel : $goodMatchRight;
		}

		#if flag ==2, good quality

		my $currPos = 0;
		foreach my $item (@token)
		{
			my $t = $item->{'t'};
			my $f = $item->{'f'};

			next if $t eq '';
			
			if ($t =~/^\d+$/)
			{
				$currPos += $t if $t> 0;
				
			}
			elsif ($t=~/^[A|C|G|T|N]$/i)
			{
				#substitution, always single nucleotide
				my $refBase = $t;
				#my $altBase = substr ($SEQ, $currPos, 1);#currPos does not include insertions, so the altBase might not be correct
				push @{$MDitems{'X'}},{pos=>$currPos, refBase=>$refBase, flag=>$f};
				$currPos+=1;
			}
			elsif ($t=~/^\^(\S+)$/)
			{
				#deletion occurs between $currPos-1 and $currPos in the read
				my $refBase = $1;
				#my $altBase = ".";
			 	$MDitems{'D'}{$currPos} = {pos=>$currPos, refBase=>$refBase, flag=>$f};
			}
			else
			{
				Carp::croak "something is wrong in $MDStr: $t\n";
			}
		}
		#print "MD=$MDStr, items=", Dumper (\%MDitems);
	}
	return \%MDitems;
}


#return 0 if no alignment
sub samToBedAndMutation
{
	my $sam = $_[0];
	return 0 if $sam->{"CIGAR"} eq '*'; #no alignment
	
	my $flagInfo = decodeSamFlag ($sam->{"FLAG"});
	
	my $strand = $flagInfo->{'query_strand'};
	
	my $name = $sam->{"QNAME"};
	my $chrom = $sam->{"RNAME"};
	my $chromStart = $sam->{"POS"} - 1;
	my $CIGAR = $sam->{"CIGAR"};
	my $SEQ = $sam->{"SEQ"};

	my $TAGS = $sam->{"TAGS"} ? $sam->{"TAGS"} : "";
	

	my $score = 0;
	if ($TAGS=~/NM\:\S*\:(\d+)/)
	{
		$score = $1;
	}

	my $matchStart = 1;

	#remove hard cliped nucleotides
	if ($CIGAR =~/^\d+H(.*?)$/)
	{
		$CIGAR = $1;
		$matchStart = 0 if $strand eq '+';
		
	}
	
	if ($CIGAR =~/^(.*?)\d+H$/)
	{
		$CIGAR = $1;
		$matchStart = 0 if $strand eq '-';
	}
	
	#remove soft cliped nucleotides
	if ($CIGAR =~/^(\d+)S(.*?)$/)
	{
		my $size = $1;
		$CIGAR = $2;
		$SEQ = substr ($SEQ, $size); #trim left
		$matchStart = 0 if $strand eq '+';
	}

	if ($CIGAR =~/^(.*?)(\d+)S$/)
	{
		$CIGAR = $1;
		my $size = $2;
		$SEQ = substr ($SEQ, 0, length($SEQ) - $size);
		$matchStart = 0 if $strand eq '-';
	}

	#print join ("\t", $name, $CIGAR), "\n";
	#deal with the rest
	if ($CIGAR=~/[^\d+|M|N|I|D]/g)
	{
		Carp::croak "unexpected CIGAR string: $CIGAR in $name: $SEQ\n";
	}

	my $MDitems = parseMD ($TAGS, $CIGAR, $SEQ);


	my (@blockSizes, @blockStarts);
	
	my $currPosRef = 0; #tracks the current genomic position (relative)
	my $currPosRead = 0; #tracks the current read position

	my $extendBlock = 0;
	my @mutations;

	my $currMismatchIdx = 0;
	
	#note in sam, the sequence reported is always the sequence on the positive strand of the chromosome, not necessarily the read sequence
	#the coordinates are also relative to the positive strand
	my $inDel = 0;
	my $ins = 0; #track the number of inserted nucleotides in cigar string, note that these are not counted in the MD tags	

	while ($CIGAR =~/(\d+)([M|N|I|D])/g)
	{
		my ($size, $type) = ($1, $2);
		if ($type eq 'I' || $type eq 'D')
		{
			$inDel += $size;
			$extendBlock = 1;
			
			if ($type eq 'I')
			{
				my $mutationStart = $chromStart + $currPosRef;
				#insertion in read occurs between $currPosRef - 1 and $currPosRef

				my $mutationEnd = $mutationStart;
				my $mutationShift = $currPosRead;
				
				#very difficult to track for insertions, so we just use the position of the insertion in the read as an approximation
				my $flag = 0;
				$flag = 2 if $mutationShift >= $minAnchorIndel && $mutationShift + $size - 1 <= length($SEQ) - 1 - $minAnchorIndel;
				my $altBase = substr ($SEQ, $currPosRead, $size);

				push @mutations, {
					type=>$type,
					refBase=> '.',
					altBase=> $altBase,
					flag=>$flag,
					chrom=>$chrom,
					chromStart=>$mutationStart,
					chromEnd => $mutationEnd,
					name=>$name,
					score=>$mutationShift,
					strand=>$strand } if $altBase=~/^[ACGT]*$/i;
				$ins += $size;
				#$currPosRead += $size; #we do not shift the read position to be consistent with the MD tag
	
			}
			else
			{
				#deletion in reads
				my $mutationStart = $chromStart + $currPosRef;
				my $mutationEnd = $mutationStart + $size - 1;
				my $mutationShift = $currPosRead;
				#the deletion is between $currPosRead -1 and $currPosread

				my $refBase = $MDitems->{'D'}{$currPosRead}{'refBase'};			
				push @mutations, {
                    type=>$type,
                    refBase=> $refBase,
                    altBase=> '.',
					flag=> $MDitems->{'D'}{$currPosRead}{'flag'},
                    chrom=>$chrom,
                    chromStart=>$mutationStart,
                    chromEnd => $mutationEnd,
                    name=>$name,
                    score=>$mutationShift,
                    strand=>$strand } if $refBase =~/^[ACGT]*$/i;				

				#deal with genomic coordinates
				my $n = @blockSizes;
				if ($n < 1)
				{
					$chromStart += $size;	
				}
				else
				{
					$blockSizes[$#blockSizes] += $size;
					$currPosRef += $size;
				}
			}
		}
		elsif ($type eq 'M')
		{
			
			#substitutions have to be in the match blocks
			my $i = 0;
			for ($i = $currMismatchIdx; $i < @{$MDitems->{'X'}}; $i++)
			{

				my $m = $MDitems->{'X'}[$i];
				#print "currPosRead = $currPosRead, size=$size\n", Dumper ($m), "\n";

				if ($m->{'pos'} < $currPosRead + $size)
				{
					#it is in this block
					my $mutationStart = $chromStart + $currPosRef + $m->{'pos'} - $currPosRead;
					my $altBase = substr ($SEQ, $m->{'pos'} + $ins, 1);
					push @mutations, {
						type=>'X',
						refBase=>$m->{'refBase'},
						altBase=>$altBase,
						#altBase=>$m->{'altBase'},
						flag=>$m->{'flag'},
						chrom=>$chrom,
						chromStart=>$mutationStart,
						chromEnd=>$mutationStart,
						name=>$name,
						score=>$m->{'pos'},
						strand=>$strand} if $m->{'refBase'}=~/^[ACTG]*$/i && $altBase=~/^[ACTG]*$/i;
				}
				else
				{
					last;
				}
			}
			$currMismatchIdx = $i;

			#deal with genomic coordinates
			if ($extendBlock && @blockSizes > 0)
			{
				#extend the previous block
				my $n = @blockSizes;
				$blockSizes[$n-1] += $size;
			}
			else
			{
				push @blockSizes, $size;
				push @blockStarts, $currPosRef;
			}
			$extendBlock = 0;

			$currPosRef += $size;
			$currPosRead += $size;
		}
		elsif ($type eq 'N')
		{
			$currPosRef += $size;
		}
		else
		{
			Carp::croak "wrong type $type\n";
		}
	}

	my $blockCount = @blockSizes;
	my $chromEnd = $chromStart + $blockStarts[$blockCount-1] + $blockSizes[$blockCount-1] - 1;
	$score -= $inDel unless $inDelInScore;	
	Carp::croak "negative score:", Dumper ($sam), "\n" unless $score >= 0;
	#we count substitutions as mismatches (caused by sequencing errors) assuming indel are caused by crosslinking

	my $bed = {
		chrom=>$chrom,
		chromStart=>$chromStart,
		chromEnd=>$chromEnd,
		name=>$name,
		score=>$score,
		strand=>$strand,
		thickStart=>$chromStart,
		thickEnd=>$chromEnd,
		itemRgb=>0,
		blockCount=>$blockCount,
		blockSizes=>\@blockSizes,
		blockStarts=>\@blockStarts,
		matchStart=>$matchStart
	};

	return {tag=>$bed, mutations=>\@mutations};
}	

