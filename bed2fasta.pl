#!/usr/bin/perl -w
#

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use Carp;
use Cwd;
use Bio::SeqIO;
use Data::Dumper;

use MyConfig;
use Common;
use Bed;
use Sequence;

my $prog = basename($0);

my $maskRepeats = 0;
my $neighbor = ''; # or 'up' or 'down'
my ($left, $right) = (0, 0);

my $twoBitFile = "";

my $proj = "/home5/zhang1/zhangc/proj";
my $twoBitToFa = "twoBitToFa";
my $nibFrag = "nibFrag";
my $RepeatMasker = "RepeatMasker"; #"/info/bin/repeatmasker/RepeatMasker";
my $organism = "mm9";
my $genomeDir = ""; #"/data/zhangc/data/genomes/mm6";
my $cache = getDefaultCache($prog);
my $keepCache = 0;
my $singleLine = 0;
my $chrLen = {};
my $chrLenFile = "";
my $capital = 0;
my $mutate = 0;
my $verbose = 0;

my $outBedFile = "";

my @ARGV0 = @ARGV;

GetOptions (
			'n|neighbor:s'=>\$neighbor,
			'l|left:i'=>\$left,
			'r|right:i'=>\$right,
			'ob:s'=>\$outBedFile,	
			'c|cache:s'=>\$cache,
			'keep-cache'=>\$keepCache,
			'org:s'=>\$organism,
			'gd:s'=>\$genomeDir,
			'twoBitFile:s'=>\$twoBitFile,
			'chrLen:s'=>\$chrLenFile,
			'nibFrag:s'=>\$nibFrag,
			'twoBitToFa:s'=>\$twoBitToFa,
			'mutate'=>\$mutate,
			'capital'=>\$capital,
			's|single-line'=>\$singleLine,
#			'm|mask_repeat'=>\$maskRepeats,
#			'rm:s'=>\$RepeatMasker,
			'v|verbose'=>\$verbose
);

if (@ARGV != 2)
{
	print "Extract sequences specified by a bed file\n";
	print "Usage: $prog [options] <mutation.bed> <mutation.fa>\n";
	print "OPTIONS:\n";
	print " -n  [string]: get neighbor sequences <up|down>\n";
	print " -l  [int]   : extension on the left, with sign <$left>\n";
	print " -r  [int]   : extension on the right, with sign <$right>\n";
	print " -mutate     : introduce random mutations with the number given in the score column\n";
	#print " -m          : mask repeats <on|[off]> (obsolete feature)\n";
	print " -s          : print each sequence in a single line\n";
	print " -capital    : print capital letters\n";
	print " -ob: output bed file name, effective only with -n\n";
	print " -v : verbose\n";
	print "External programs and dir\n";
	print " -nibFrag    [string]: path to the nibFrag program [$nibFrag]\n";
	print " -twoBitToFa [string]: path to the twoBitToFa program [$twoBitToFa]\n";
	#print " -rm         [string]: path to the RepeatMasker program [$RepeatMasker]\n";
	print " -c          [string]: cache dir [$cache]\n";
	print " --keep-cache        : keep cache when the job is done\n";
	print " -org        [string]: organism <[mm9]|mm8|mm6|mm5|rn4|rn3|hg17|hg18|sacCer1|dm2|ce2>\n";
	print " -gd         [string]: directory of nib files [$genomeDir] (will override -org)\n";
	print " -twoBitFile [string]: twoBitFile [$twoBitFile] (will override -org or -gd)\n";
	print " -chrLen     [string]: chromosome length file\n";
	print "\n";
	exit (1);
}


print "CMD=$prog ", join(" ", @ARGV0), "\n" if $verbose;

my ($in, $out) = @ARGV;

$cache = getFullPath ($cache);


#if genomeDir is specified, use that
#

if ($twoBitFile ne '')
{
	Carp::croak "$twoBitFile does not exists\n" unless -f $twoBitFile;
	$twoBitFile = getFullPath ($twoBitFile);
}
elsif ($genomeDir ne '')
{
	Carp::croak "$genomeDir does not exists\n" unless -d $genomeDir;
	$genomeDir = getFullPath ($genomeDir);
}
elsif ($organism ne '')
{
	#if organism is specified, use that

	#Carp::croak "No species specified\n" if $organism eq '' && $twoBitFile eq '';
	$genomeDir = MyConfig::getGenomeDir($organism);
	#my $chrLenFile = "$genomeDir/chrLen.txt";
	
	$chrLenFile = "$genomeDir/chrLen.txt";
	$twoBitFile = "$genomeDir/$organism.2bit";

	$genomeDir .="/nib";
	Carp::croak "neither $genomeDir nor $twoBitFile exists\n" unless (-d $genomeDir || -f $twoBitFile);
}
else
{
	Carp::croak "No species specified\n";
}

=obsolete
sub maskRepeat
{
	my ($in, $out, $cache) = @_;
	$in = Common::getFullPath ($in);
	$out = Common::getFullPath ($out);
	$cache = Common::getFullPath ($cache);
	
	my $RMdir = "$cache/RM" .rand ();
	system ("mkdir $RMdir");

	#change directory because RM generate output in the working directory
	my $workDir = cwd;
	chdir ($RMdir);
	
	#print cwd, "\n";
	#print "rm cmd = $RepeatMasker -dir $RMdir $fastaTmp &>/dev/null\n";
	
	my $rmCmd = "$RepeatMasker -dir $RMdir $in &>/dev/null";
	print "$prog: $rmCmd\n" if $verbose;
	system ($rmCmd);
	chdir ($workDir);
	
	my $fastaTmpRm = "$RMdir/" .basename ($in) .".masked";
	system ("mv $fastaTmpRm $out");
	system ("rm -rf $RMdir");
}
=cut



$outBedFile = '' if $outBedFile eq 'None';

#if (($left || $right) && $outBedFile eq '')
#{
#	$outBedFile = "$in.shift";
#}


print "loading chromosome length from $chrLenFile\n" if (-f $chrLenFile) && $verbose;

if ($chrLenFile)
{
	my $fin;
	open ($fin, "<$chrLenFile") || Carp::croak "can not open file $chrLenFile to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
			
		my ($chrom, $len) = split ("\t", $line);
		$chrLen->{$chrom} = $len;
	}
	close ($fin);
}


print "loading regions from $in ...\n" if $verbose;
my $regions = readBedFile($in, $verbose);

my $rmCache = 0;
if (not (-d $cache))
{
	system ("mkdir $cache");
	$rmCache = 1;
}

print "extracting sequences to $out ...\n" if $verbose;

#my $fout = new FileHandle;
#open ($fout, ">$out") || Carp::croak "can not open file $out to write\n";
my $seqOutStream = Bio::SeqIO->new (-file=> ">$out", -format=>'fasta');


my $foutBed;
if ($outBedFile ne '')
{
	#Carp::croak "$outBedFile already exists\n" if -f $outBedFile;
	open ($foutBed, ">$outBedFile") || Carp::croak "can not open file $outBedFile to write\n";
}

for (my $i = 0; $i < @$regions; $i++)
{

	print "$i ...\n" if $verbose && $i % 5000 == 0;	
	my $chrom = $regions->[$i]->{"chrom"};
	my $chromStart = $regions->[$i]->{"chromStart"};
	my $chromEnd = $regions->[$i]->{"chromEnd"};
	my $strand = '+';
	if (exists $regions->[$i]->{"strand"})
	{
		$strand = '-' if $regions->[$i]->{"strand"} eq '-';
	}

	#modified to allow extension
	if ($neighbor eq '')
	{
		if ($strand eq '+')
		{
			$chromStart += $left;
			$chromEnd += $right;
		}
		else
		{
			$chromStart -= $right;
			$chromEnd -= $left;
		}
		if ($chromEnd < $chromStart)
		{
			print "$chrom:$chromStart-$chromEnd is too small for truncation\n" if $verbose;
			next;
		}
	}
	elsif ($neighbor eq 'down' && $right - $left > 0)
	{
		if ($strand eq '+')
		{
			$chromStart = $chromEnd + $left;
			$chromEnd = $chromEnd + $right;
		}
		else
		{
			$chromEnd = $chromStart - $left;
			$chromStart = $chromStart - $right;
		}
	}
	elsif ($neighbor eq 'up' && $right - $left > 0)
	{
		if ($strand eq '+')
		{
			$chromEnd = $chromStart + $right;
			$chromStart = $chromStart + $left;
		}
		else
		{
			$chromStart = $chromEnd - $right;
			$chromEnd = $chromEnd - $left;
		}
	}
	else
	{
		Carp::croak "which region do you want?\n";
	}

	$chromEnd += 1; #to consistent with nibFrag and bed file format
	
	my ($upt, $dnt) = (0, 0); #if the left or right end reach the end of the chromosome
	if ($chromStart < 0)
	{
		#number of nucleotides truncated upstream
		if ($strand eq '+')
		{
			$upt = -$chromStart;
		}
		else
		{
			$dnt = -$chromStart;
		}
		$chromStart = 0;
	}
	
	if ($chromEnd < 1)
	{
		warn "negative chromEnd: ", Dumper ($regions->[$i]), "\n";
		next;
		#$chromEnd = 1;
	}

	if (ref ($chrLen) && exists $chrLen->{$chrom} && $chromEnd > $chrLen->{$chrom})
	{
		#number of nucleotides truncated downstream
		if ($strand eq '+')
		{
			$dnt = $chromEnd - $chrLen->{$chrom};
		}
		else
		{
			$upt = $chromEnd - $chrLen->{$chrom};
		}
		$chromEnd = $chrLen->{$chrom};
		#$dnt = 1;
	}

	if (ref ($chrLen) && exists $chrLen->{$chrom} && $chromStart >= $chrLen->{$chrom})
	{
		warn "chromStart goes beyond the end of the chromosome:", Dumper ($regions->[$i]), "\n";
		next;
		#$chromStart = $chrLen->{$chrom} - 1;
	}

	my $nibFragTmp = $cache . "/" . join ("_", $chrom, $chromStart, $chromEnd) . ".fa." . rand();
	my $header = "$chrom:$chromStart-$chromEnd";
	if (exists $regions->[$i]->{"name"})
	{
		if (length ($regions->[$i]->{"name"})> 0)
		{
			$header = $regions->[$i]->{"name"};
			$header .= ":$neighbor" if ($neighbor ne '');	#to generate unique names
		}
	}
	my $id = $header;
	my $name = $header;
	
	if ($upt || $dnt)
	{
		print "$name reach the end of the chromosome\n" if $verbose;
	}

	my $desc = "/strand=$strand";
	$desc .= "\t/upt=$upt" if $upt;
	$desc .= "\t/dnt=$dnt" if $dnt;

	$header .= "\t$desc";


	if (-f $twoBitFile)
	{
		#NOTE: the strand is not reversed here, need to do this later
		my $cmd = "$twoBitToFa -seq=$chrom -start=$chromStart -end=$chromEnd $twoBitFile $nibFragTmp";
		my $ret = system ($cmd);
		Carp::croak "Command [$cmd] crashed: $?\n" unless $ret == 0;
	}
	elsif (-d $genomeDir)
	{
		#get sequence by nibFrag
		$strand = 'm' if $strand eq '-';
	
		my $chromNib = "$genomeDir/$chrom.fa.nib";
		$chromNib = "$genomeDir/$chrom.nib" unless (-f $chromNib);

		Carp::croak "$chromNib does not exists\n" unless -f $chromNib;
	
		my $cmd = "$nibFrag -masked -name=\"$header\" $chromNib $chromStart $chromEnd $strand $nibFragTmp";
	
		#print $cmd, "\n";
		my $ret = system ($cmd);
		Carp::croak "Command [$cmd] crashed: $?\n" unless $ret == 0;
		
		#my $fin = new FileHandle;
		#open ($fin, "<$nibFragTmp") || Carp::croak "can not open file $nibFragTmp to read\n";
		#my @lines = <$fin>;
		#close ($fin);
		#unlink $nibFragTmp;
		#print $fout @lines;
	}
	else
	{
		Carp::croak "Something wrong\n";	
	}

	
	#my $fin = new FileHandle;
	#open ($fin, "<$nibFragTmp") || Carp::croak "can not open file $nibFragTmp to read\n";
	#my @lines;
	#while (my $line = <$fin>)
	#{
	#	chomp $line;
	#	next if $line =~/^\s*$/;
	#	next if $line=~/^\>/;
	#	push @lines, $line;
	#}
	#close ($fin);
	#
	
	my $seqInStream = Bio::SeqIO->new (-file => $nibFragTmp, -format=> 'fasta');
	my $seq = $seqInStream->next_seq ();

	$seq->id ($id);
	$seq->desc ($desc);
	if ((-f $twoBitFile) && $strand eq '-')
	{
		#now we reverse the strand
		$seq = $seq->revcom ();
	}

	unlink $nibFragTmp;
	#print $fout ">$header\n";
	#print $fout join ("\n", @lines), "\n";
	
	if ($mutate)
	{
		#introduce random mutations
		Carp::croak "Number of mutations must be specified in the score column\n" unless exists $regions->[$i]->{"score"};
		my $numMutation = $regions->[$i]->{"score"};
		
		if ($numMutation > 0)
		{
			my $mutatedSeq = mutateSeq ($seq->seq(), $numMutation);
			$seq->seq($mutatedSeq);
		}
	}

	if ($capital)
	{
		my $seqStr = $seq->seq ();
		$seqStr =~tr/a-z/A-Z/;
		$seq->seq($seqStr);
	}

	$seqOutStream->write_seq ($seq);

	if (-f $outBedFile) #the file should have been created at this point if necessary
	{
		my $r = {chrom=>$regions->[$i]->{"chrom"},
				chromStart=> $chromStart,
				chromEnd=> $chromEnd-1,
				name=>$name,
				strand=> ($strand eq '+')? '+' : '-'
		};
		$r->{"score"} = (exists $regions->[$i]->{"score"})? $regions->[$i]->{"score"} : 0;

		print $foutBed bedToLine ($r), "\n";
	}
}
#close ($fout);

if ($outBedFile ne '')
{
	close ($foutBed);
}

=obsolete
if ($maskRepeats)
{
	print "masking repeats\n" if $verbose;
	maskRepeat ($out, $out, $cache);	
}
=cut

if ($singleLine)
{
	print "convert to single line format ...\n" if $verbose;
	my $outSingle = "$out.single";
	my $fout;
	open ($fout, ">$outSingle") || Carp::croak "can not open file $outSingle to write\n";
	my $seqIO = Bio::SeqIO->new (-file=>$out, -format=>'Fasta');
	while (my $seq = $seqIO->next_seq())
	{
		my $header = $seq->id();
		if ($seq->desc() ne '')
		{
			$header .= "\t" . $seq->desc();
		}
		print $fout ">$header\n";
		print $fout $seq->seq(), "\n";
	}
	close ($fout);
	system ("mv $outSingle $out");
}

system ("rm -rf $cache") unless $keepCache;


