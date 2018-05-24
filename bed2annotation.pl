#!/usr/bin/perl -w


use strict;
use Getopt::Long;

use Carp;
use File::Basename;
use Data::Dumper;

use MyConfig;

my $prog = basename ($0);
my $cmdDir = dirname ($0);
my $verbose = 0;

my @ARGV0 = @ARGV;


my $big = 0;
my $confFile = "";
my $dbkey = "";
my $cache = getDefaultCache ($prog);
my $keepCache = 0;

my $separateStrand = 0;

my $get_gene = 0;
my $get_rmsk = 0;
my $get_miRNA = 0;
my $get_region = 0;
my $get_custom = ""; #custom features, a bed file
my $customName = "custom_feature";
my $customSummaryMethod = "all";

my $summaryFile = "";

GetOptions (
	'conf=s'=>\$confFile,
	'ss'=>\$separateStrand,
	'dbkey=s'=>\$dbkey,
	'big'=>\$big,
	'gene'=>\$get_gene,
	'rmsk'=>\$get_rmsk,
	'miRNA'=>\$get_miRNA,
	'region'=>\$get_region,
	'custom:s'=>\$get_custom,
	'custom-name:s'=>\$customName,
	'custom-summary:s'=>\$customSummaryMethod,
	'c:s'=>\$cache,
	'keep-cache'=>\$keepCache,
	'summary:s'=>\$summaryFile,
	'v'=>\$verbose);


if (@ARGV != 2)
{
	print "add various annotations to bed regions\n";

	print "Usage: $prog [options] <in.bed> <out.txt>\n";
	print " -conf            [file]   : configuration file with input datasets\n";
	print " -dbkey           [string] : genome build name (hg19|mm10)\n";
	print " -ss                       : consider the two strands separately when possible\n";
	print " -big                      : big file\n";
	print " -gene                     : annotate overlapping gene (id and symbol)\n";
	print " -rmsk                     : annotate overlapping RepeatMasked sequences (type and \%)\n";
	print " -miRNA                    : annotate miRNA (miRNA_name)\n";
	print " -region                   : annotate genomic breakdown\n";
	print " -custom          [file]   : annotate custom features in the provided BED file\n";
	print " --custom-name    [string] : naming the custom feature\n";
	print " --custom-summary [string] : method to summarize custom annotation ([all]|max_num|min_num|max_overlap)\n";
	print " -summary         [file]   : print summary information\n";
	print " -c               [string] : cache dir ($cache) \n";
	print " --keep-cache              : keep cache when the job is done\n";
	print " -v                        : verbose\n";
	exit (1);
}

my ($inBedFile, $outFile) = @ARGV;
my $msgio = $outFile eq '-' ? *STDERR :  *STDOUT;
print $msgio "CMD = $prog ", join (' ', @ARGV0), "\n" if $verbose;

system ("mkdir $cache");

my %analyses;

my $verboseFlag = $verbose ? '-v' : '';
my $bigFlag = $big ? '-big' : '';
my $keepCacheFlag = $keepCache ? '--keep-cache' : '';


Carp::croak "no annotations chosen\n" if $get_gene + $get_rmsk + $get_miRNA + $get_region == 0 && $get_custom eq '';

$analyses{'gene'} = 1 if $get_gene;
$analyses{'rmsk'} = 1 if $get_rmsk;
$analyses{'miRNA'} = 1 if $get_miRNA;
$analyses{'region'} = 1 if $get_region;
$analyses{'custom'} = 1 if $get_custom ne '';


if ($confFile ne '')
{
	Carp::croak "$confFile does not exist\n" unless -f $confFile;
}
else
{
	#take the default, assuming it is located in the same directory as the script
	$confFile = "$cmdDir/ctk.loc";
}


my $customBedFile = $get_custom;

if ($customBedFile ne '')
{
	Carp::croak "cannot file $customBedFile" unless -f $customBedFile;
}


#save the original id
my $inBedFileTmp = "$cache/in.bed";
my $inIdFile = "$cache/in.id";
my $lineNoFile = "$cache/in.lineNo";
#Carp::croak "OK\n";

#NOTE: "cat: write error: Broken pipe" means that cat was attempting to write stdout to the stdin of another process but the other process terminated. This will happen for every one of the cat <some file> | head -<n> commands when there are more than <n> lines in <some file>. The crude solution would be to discard the error messages

#my $cmd = "cat $inBedFile | grep -v \"^track\" 2>/dev/null | head -n 1 | awk '{print NF}' > $lineNoFile";

my $cmd = "cat $inBedFile 2>/dev/null | grep -v \"^track\" 2>/dev/null | head -n 1 | awk '{print NF}' > $lineNoFile";

system ($cmd);

my $colNum = `cat $lineNoFile`;
chomp $colNum;

Carp::croak "The input Bed file must have at least 3 columns\n" unless $colNum >= 3;

if ($colNum ==3)
{
	#no id, so use numbers
	my $cmd = "grep -v \"^track\" $inBedFile | awk 'BEGIN{i=0} {i=i+1; print i}' > $inIdFile";
	system ($cmd);
}
else
{
	#keep the id
	my $cmd = "grep -v  \"^track\" $inBedFile | cut -f 4 > $inIdFile"; 
	system ($cmd);
}

#replace id in case there is duplicates in the original id
$cmd = "grep -v \"^track\" $inBedFile | awk 'BEGIN{i=0} {i=i+1; if(NF<6) {print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"i} else {print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"i\"\\t\"\$5\"\\t\"\$6}}' > $inBedFileTmp";
print $cmd, "\n" if $verbose;
system ($cmd);

#no header now
$inBedFile = $inBedFileTmp;
my $totalN = `wc -l $inBedFile | awk '{print \$1}'`; chomp $totalN;

my $fout_summary;
open ($fout_summary, ">$summaryFile") if $summaryFile ne '';

my $date = `date`; chomp $date;
my $host = `hostname`; chomp $host;

print $fout_summary "$date at $host\n\n" if $summaryFile ne '';
foreach my $analysis (sort keys %analyses)
{
	print "get $analysis annotation ...\n" if $verbose;
	my $locationInfo = $analysis eq 'custom' ? {custom=>$customBedFile} : getLocationInfo ($confFile, $dbkey, $analysis);

	if (keys %$locationInfo == 0)
	{
		my $msg = "cannot locate $analysis annotation bed file";
		$msg .= " for $dbkey" if $dbkey ne '';

		Carp::croak $msg, "\n";
	}

	#Carp:croak Dumper ($locationInfo), "\n";
	my $tmpOutFile = "$cache/$analysis.out";

	if (-f $summaryFile)
	{
		print $fout_summary "\n\nSummary of $analysis annotation\n";
		print $fout_summary "-"x30, "\n\n";
	}

	if ($analysis eq 'gene')
	{
		my $transcriptBedFile = $locationInfo->{'transcript'};
		my $transcript2geneFile = $locationInfo->{'transcript2gene'};
		my $gene2symbolFile = $locationInfo->{'gene2symbol'};

		my $ssFlag = $separateStrand ? '-ss' : '';

		#gene_id gene_symbol
		my $cmd = "perl $cmdDir/bed2gene2symbol.pl --no-region-id $bigFlag $ssFlag $verboseFlag $inBedFile $transcriptBedFile $transcript2geneFile $gene2symbolFile $tmpOutFile";
		print $msgio "$cmd\n" if $verbose;

		my $ret = system ($cmd);
		Carp::croak "$cmd failed: $?\n" if $ret != 0;

		if (-f $summaryFile)
		{
			my $geneN = `awk '{if(NF>1) {print \$0}}' $tmpOutFile | wc -l | awk '{print \$1}'`; chomp $geneN;
			#Carp::croak"geneN = $geneN\n";

			print $fout_summary join("\t", "Region", "No.", "%"), "\n";
			print $fout_summary join ("\t", "Intervals overlapping with UCSC/RefSeq genes", $geneN, sprintf ("%.1f", $geneN / $totalN * 100)), "\n";
			print $fout_summary join ("\t", "Other", $totalN - $geneN, sprintf ("%.1f", 100 - $geneN / $totalN * 100)), "\n";
		}
	}
	elsif ($analysis eq 'rmsk')
	{
		my $rmskBedFile = $locationInfo->{'rmsk'};
		my $ssFlag = "";

		bed2annot ($inBedFile, $rmskBedFile, $tmpOutFile, "max_overlap", $analysis, $cache, $ssFlag, $bigFlag, $verboseFlag, $keepCacheFlag);
		if (-f $summaryFile)
		{
			my $rmskN = `awk '{if(NF>1) {print \$0}}' $tmpOutFile | wc -l | awk '{print \$1}'`; chomp $rmskN;
			print $fout_summary join("\t", "Region", "No.", "%"), "\n";
			print $fout_summary join ("\t", "Intervals overlapping repeat masked region", $rmskN, sprintf ("%.1f", $rmskN / $totalN * 100)), "\n";
			print $fout_summary join ("\t", "Other", $totalN - $rmskN, sprintf ("%.1f", 100 - $rmskN / $totalN * 100)), "\n";
		}

	}
	elsif ($analysis eq 'miRNA' || $analysis eq 'custom')
	{
		my $summaryMethod = $analysis eq 'custom' ? $customSummaryMethod : 'all';

		my $featureBedFile = $locationInfo->{$analysis};
		my $ssFlag = $separateStrand ? '-ss' : '';
	
		bed2annot ($inBedFile, $featureBedFile, $tmpOutFile, $summaryMethod, $analysis, $cache, $ssFlag, $bigFlag, $verboseFlag, $keepCacheFlag);

		if (-f $summaryFile)
		{
			my $featureName = $analysis eq 'custom' ? $customName : $analysis;
			my $featureN = `awk '{if(NF>0) {print \$0}}' $tmpOutFile | wc -l | awk '{print \$1}'`; chomp $featureN;
			print $fout_summary join("\t", "Region", "No.", "%"), "\n";
			print $fout_summary join ("\t", "Intervals overlapping $featureName", $featureN, sprintf ("%.1f", $featureN / $totalN * 100)), "\n";
			print $fout_summary join ("\t", "Other", $totalN - $featureN, sprintf ("%.1f", 100 - $featureN / $totalN * 100)), "\n";
		}
	}
	elsif ($analysis eq 'region')
	{
		
		my $cache2 = "$cache/region_cache";
		system ("mkdir $cache2");

		my $genicBedFile = $locationInfo->{'genic'};
		my $upstreamBedFile = $locationInfo->{'upstream_10k'};
		my $downstreamBedFile = $locationInfo->{'downstream_10k'};
		my $exonBedFile = $locationInfo->{'exon'};
		my $exon_5utrBedFile = $locationInfo->{'5utr'};
		my $exon_3utrBedFile = $locationInfo->{'3utr'};
		my $exon_cdsBedFile = $locationInfo->{'cds'};

		my $ssFlag = $separateStrand ? '-ss' : '';

		my $bed_vs_exonBedFile = "$cache2/bed_vs_exon.bed";
		my $bed_vs_exon_5utrBedFile = "$cache2/bed_vs_5utr_exon.bed";
		my $bed_vs_exon_3utrBedFile = "$cache2/bed_vs_3utr_exon.bed";
		my $bed_vs_exon_cdsBedFile = "$cache2/bed_vs_cds_exon.bed";

		my $bed_vs_genicBedFile = "$cache2/bed_vs_genic.bed";
		my $bed_vs_upstreamBedFile = "$cache2/bed_vs_upstream.bed";
		my $bed_vs_downstreamBedFile = "$cache2/bed_vs_downstream.bed";

		my @fileInfo = (
			[$genicBedFile, $bed_vs_genicBedFile],
			[$upstreamBedFile, $bed_vs_upstreamBedFile],
			[$downstreamBedFile, $bed_vs_downstreamBedFile],
			[$exonBedFile, $bed_vs_exonBedFile],
			[$exon_5utrBedFile, $bed_vs_exon_5utrBedFile],
			[$exon_3utrBedFile, $bed_vs_exon_3utrBedFile],
			[$exon_cdsBedFile, $bed_vs_exon_cdsBedFile]);

		foreach my $g (@fileInfo)
		{
			my $in = $g->[0];
			my $out = $g->[1];
			my $cmd = "perl $cmdDir/tagoverlap.pl $bigFlag $ssFlag $verboseFlag $keepCacheFlag -region $in $inBedFile $out";
			my $ret = system ($cmd);
			Carp::croak "$cmd failed: $?\n" if $ret != 0;
		}

		my $bed_vs_genic_idFile = "$cache2/bed_vs_genic.id";
		my $bed_vs_upstream_idFile = "$cache2/bed_vs_upstream.id";
		my $bed_vs_downstream_idFile = "$cache2/bed_vs_downstream.id";
		my $bed_vs_exon_idFile = "$cache2/bed_vs_exon.id";

		@fileInfo = (
			[$bed_vs_genicBedFile, $bed_vs_genic_idFile],
			[$bed_vs_upstreamBedFile, $bed_vs_upstream_idFile],
			[$bed_vs_downstreamBedFile, $bed_vs_downstream_idFile],
			[$bed_vs_exonBedFile, $bed_vs_exon_idFile]);

		foreach my $g (@fileInfo)
		{
			my $in = $g->[0];
			my $out = $g->[1];
			my $cmd = "awk '{print \$4}' $in | awk -F \"//\" '{print \$1}' | sort | uniq > $out";
			my $ret = system ($cmd);
			Carp::croak "$cmd failed: $?\n" if $ret != 0;
		}

		my $bed_vs_genic_and_upstream_idFile = "$cache2/bed_vs_genic_and_upstream.id";
		my $bed_vs_genic_and_downstream_idFile = "$cache2/bed_vs_genic_and_downstream.id";

		my $cmd = "join $bed_vs_genic_idFile  $bed_vs_upstream_idFile > $bed_vs_genic_and_upstream_idFile";
		system ($cmd);

		$cmd = "join $bed_vs_genic_idFile  $bed_vs_downstream_idFile > $bed_vs_genic_and_downstream_idFile";
		system ($cmd);

		my $bed_vs_exon_5utr_idpairFile = "$cache2/bed_vs_5utr_exon.idpair";
		my $bed_vs_exon_3utr_idpairFile = "$cache2/bed_vs_3utr_exon.idpair";
		my $bed_vs_exon_cds_idpairFile = "$cache2/bed_vs_cds_exon.idpair";

		@fileInfo = (
			[$bed_vs_exon_5utrBedFile, $bed_vs_exon_5utr_idpairFile, "5'\\''UTR"],
			[$bed_vs_exon_3utrBedFile, $bed_vs_exon_3utr_idpairFile, "3'\\''UTR"],
			[$bed_vs_exon_cdsBedFile, $bed_vs_exon_cds_idpairFile, "CDS"]);
		
		foreach my $g (@fileInfo)
		{
			my $in = $g->[0];
			my $out = $g->[1];
			my $type = $g->[2];
			my $cmd = "awk '{print \$4}' $in | awk -F \"//\" '{print \$1}' | sort | uniq | awk '{print \$1\"\\t$type\"}' > $out";
			my $ret = system ($cmd);
			Carp::croak "$cmd failed: $?\n" if $ret != 0;
		}

		my $bed_vs_intron_idpairFile = "$cache2/bed_vs_intron.idpair";
		$cmd = "perl $cmdDir/removeRow.pl $bed_vs_genic_idFile $bed_vs_exon_idFile | awk '{print \$1\"\\tintron\"}' > $bed_vs_intron_idpairFile";
		my $ret = system ($cmd);
		Carp::croak "$cmd failed: $?\n" if $ret != 0;

		my $bed_vs_exon_annotated_idFile = "$cache2/bed_vs_exon_annotated.id";

		$cmd = "cat $bed_vs_exon_5utr_idpairFile $bed_vs_exon_3utr_idpairFile $bed_vs_exon_cds_idpairFile | sort | uniq | awk '{print \$1}' > $bed_vs_exon_annotated_idFile";
		$ret = system ($cmd);
		Carp::croak "$cmd failed: $?\n" if $ret != 0;

		my $bed_vs_exon_unclassified_idpairFile = "$cache2/bed_vs_exon_unclassified.idpair";
		my $bed_vs_downstream_intergenic_idpairFile = "$cache2/bed_vs_downstream_intergenic.idpair";
		my $bed_vs_upstream_intergenic_idpairFile = "$cache2/bed_vs_upstream_intergenic.idpair";


		@fileInfo = (
			[$bed_vs_exon_idFile, $bed_vs_exon_annotated_idFile, $bed_vs_exon_unclassified_idpairFile, "exon_unclassified"],
			[$bed_vs_downstream_idFile, $bed_vs_genic_idFile, $bed_vs_downstream_intergenic_idpairFile, "downstream_10k"],
			[$bed_vs_upstream_idFile,   $bed_vs_genic_idFile, $bed_vs_upstream_intergenic_idpairFile,   "upstream_10k"]);
		
		for my $g (@fileInfo)
		{
			my $in1 = $g->[0];
			my $in2 = $g->[1];
			my $out = $g->[2];
			my $type = $g->[3];
			
			my $cmd = "perl $cmdDir/removeRow.pl $in1 $in2 | awk '{print \$1\"\\t$type\"}' > $out";
			my $ret = system ($cmd);
			Carp::croak "$cmd failed: $?\n" if $ret != 0;
		}
		
		my $bed_vs_annot_idpairFile = "$cache2/bed_vs_annot.idpair";

		$cmd = "cat $bed_vs_exon_cds_idpairFile $bed_vs_exon_5utr_idpairFile $bed_vs_exon_3utr_idpairFile $bed_vs_exon_unclassified_idpairFile ";
	  	$cmd .= "$bed_vs_intron_idpairFile $bed_vs_downstream_intergenic_idpairFile $bed_vs_upstream_intergenic_idpairFile  > $bed_vs_annot_idpairFile";
		system ($cmd);
		
		#print breakdown
		$cmd = "perl $cmdDir/selectRow.pl -f 3 -p -pt \"deep_intergenic\" -s $bed_vs_annot_idpairFile $inBedFile  | cut -f 2 > $tmpOutFile";
		$ret = system ($cmd);
		Carp::croak "$cmd failed: $?\n" if $ret != 0;

		#summary
		
		if (-f $summaryFile)
		{
			#my $totalN = `wc -l $inBedFile | awk '{print \$1}'`; chomp $totalN;
			my $genicN = `cat $bed_vs_genic_idFile $bed_vs_exon_idFile | sort | uniq | wc -l | awk '{print \$1}'`; chomp $genicN;
			
			my $exonicN = `wc -l $bed_vs_exon_idFile | awk '{print \$1}'`; chomp $exonicN;
			my $exon_cdsN = `wc -l $bed_vs_exon_cds_idpairFile | awk '{print \$1}'`; chomp $exon_cdsN;
			my $exon_5utrN = `wc -l $bed_vs_exon_5utr_idpairFile | awk '{print \$1}'`; chomp $exon_5utrN;
			my $exon_3utrN = `wc -l $bed_vs_exon_3utr_idpairFile | awk '{print \$1}'`; chomp $exon_3utrN;
			my $intronicN = $genicN - $exonicN;
			
			my $upstream_intergenicN = `wc -l  $bed_vs_upstream_intergenic_idpairFile | awk '{print \$1}'`; chomp $upstream_intergenicN;
			my $downstream_intergenicN = `wc -l $bed_vs_downstream_intergenic_idpairFile | awk '{print \$1}'`; chomp  $downstream_intergenicN;
			
			my $genic_extN = `cat  $bed_vs_genic_idFile $bed_vs_upstream_idFile $bed_vs_downstream_idFile | sort | uniq | wc -l`; chomp $genic_extN;
			my $deep_intergenicN = $totalN - $genic_extN;

			my @s = (
				[$totalN, "Total"],
				[$genicN, "Genic"],
				[$exonicN, "Exon"],
				[$exon_cdsN, "CDS exon"],
				[$exon_5utrN, "5' UTR exon"],
				[$exon_3utrN, "3' UTR exon"],
				[$intronicN, "Intron"],
				[$upstream_intergenicN, "Upstream 10K"],
				[$downstream_intergenicN, "Downstream 10K"],
				[$genic_extN, "Genic+ext10K"],
				[$deep_intergenicN, "Deep intergeic"]);

			print $fout_summary join("\t", "Region", "No.", "%"), "\n";
			foreach my $g (@s)
			{
				print $fout_summary join ("\t", $g->[1], $g->[0], sprintf ("%.1f", $g->[0] / $totalN * 100)), "\n";
			}
	
			#print "genicN=$genicN, genic_extN=$genic_extN\n";

			my $extN = $genic_extN - $genicN;
			#$upstream_intergenicN = 0;
			my $upstream_intergenicN_correct = int ($upstream_intergenicN / ($upstream_intergenicN+ $downstream_intergenicN) * $extN + 0.5) 
			if $upstream_intergenicN+ $downstream_intergenicN > 0;
			
			my $downstream_intergenicN_correct = $extN - $upstream_intergenicN_correct;


			print $fout_summary "\nPercentage of intervals in each region (adjusted):\n";

			print $fout_summary join ("\t", "Region", "%"), "\n";
			print $fout_summary join ("\t", "CDS exon", sprintf ("%.1f", $exon_cdsN * $exonicN / ($exon_5utrN + $exon_3utrN + $exon_cdsN) / $totalN * 100)), "\n";
			print $fout_summary join ("\t", "5' UTR exon", sprintf ("%.1f", $exon_5utrN * $exonicN / ($exon_5utrN + $exon_3utrN + $exon_cdsN) / $totalN * 100)), "\n";
			print $fout_summary join ("\t", "Upstream 10K", sprintf ("%.1f", $upstream_intergenicN_correct / $totalN * 100)), "\n";
			print $fout_summary join ("\t", "3' UTR exon", sprintf ("%.1f", $exon_3utrN * $exonicN / ($exon_5utrN + $exon_3utrN + $exon_cdsN) / $totalN * 100)), "\n";
			print $fout_summary join ("\t", "Downstream 10K", sprintf ("%.1f", $downstream_intergenicN_correct / $totalN * 100)), "\n";
			print $fout_summary join ("\t", "Intron", sprintf ("%.1f", $intronicN /  $totalN * 100)), "\n";
			print $fout_summary join ("\t", "Deep intergenic", sprintf ("%.1f", $deep_intergenicN /  $totalN * 100)), "\n";
		}
	}

}

close ($fout_summary) if -f $summaryFile;


#combine output

my $headerLine = "#name";
$cmd = "paste $inIdFile";

my $customColHeader = ['custom', $customName];
if ($customSummaryMethod eq 'max_overlap')
{
	$customColHeader = ['custom', $customName, $customName . "_fraction"];
}
elsif ($customSummaryMethod eq 'max_num')
{
	$customColHeader = ['custom', $customName, $customName . "_maxScore"];
}
elsif ($customSummaryMethod eq 'min_num')
{
	$customColHeader = ['custom', $customName, $customName . "_minScore"];
}
elsif ($customSummaryMethod eq 'all')
{
}
else
{
	Carp::croak "unknown summary method for custom features; $customSummaryMethod\n";
}




my @fileInfo = (
	['gene', "gene_id\tgene_symbol"],
	['rmsk', "rmsk_type\trmsk_fraction"],
	['miRNA', "miRNA"],
	['region', "region"],
	$customColHeader);

foreach my $g (@fileInfo)
{
	my @g = @$g;

	my $analysis = shift @g;
	my $head = join("\t", @g);

	next unless exists $analyses{$analysis};

	my $f = "$cache/$analysis.out";

	Carp::croak "$analysis annotation result file $f was not generate properly\n" unless -f $f;

	$headerLine .= "\t$head";
	$cmd .= " $f";
}

system ("echo \"$headerLine\" > $outFile");
$cmd .= ">> $outFile";

print $msgio $cmd, "\n" if $verbose;
system ($cmd);

system ("rm -rf $cache") unless $keepCache;
print $msgio "Done.\n" if $verbose;



sub getLocationInfo
{
	my ($conf, $dbkey, $analysis) = @_;

	my $fin;
	
	open ($fin, "<$conf") || Carp::croak "cannot open file $conf\n";
	
	my %ret;

	while (my $line = <$fin>)
	{
		chomp $line;
		next if $line =~/^\s*$/;
		next if $line =~/^\#/;

		my ($db, $ana, $path, $type) = split (/\s+/, $line);

		$path = "$cmdDir/$path" unless $path=~/^\//;
		#if a relative path is provided, we assume the annotation file is located in the same folder as the script
		#fixed by CZ, 07/31/2016

		$type = $ana unless $type;

		if ($db eq $dbkey && $ana eq $analysis)
		{
			Carp::croak "$path does not exist\n" unless -f $path;
			$ret{$type} = $path;
		}
	}
	close ($fin);
	return \%ret;
}

sub bed2annot
{
	my ($inBedFile, $annotBedFile, $outFile, $summaryMethod, $name, $cache, $ssFlag, $bigFlag, $verboseFlag, $keepCacheFlag) = @_;
	
	if ($summaryMethod eq 'max_overlap')
	{
		my $bed_vs_annot_BedFile = "$cache/bed_vs_$name.bed";

		my $cmd = "perl $cmdDir/tagoverlap.pl -d \"/##/\" $bigFlag $verboseFlag $keepCacheFlag -region $annotBedFile $inBedFile $bed_vs_annot_BedFile";
		print $msgio $cmd,"\n" if $verbose;

		my $ret = system ($cmd);
		Carp::croak "$cmd failed: $?\n" if $ret != 0;
		my $bed_vs_annot_idpairFile = "$cache/bed_vs_$name.idpair";

		$cmd = "cut -f 4,5 $bed_vs_annot_BedFile | awk -F \"/##/\" '{print \$1\"\\t\"\$2}' > $bed_vs_annot_idpairFile";
		print $msgio $cmd, "\n" if $verbose;
		system ($cmd);

		my $bed_vs_annot_idpairUniqFile = "$cache/bed_vs_$name.idpair.uniq";
		$cmd = "perl $cmdDir/uniqRow.pl $verboseFlag -value 2 -c max_num $bed_vs_annot_idpairFile $bed_vs_annot_idpairUniqFile";
		$ret = system ($cmd);
		Carp::croak "$cmd failed: $?\n" if $ret != 0;

		#print annotation_id overlap_fraction
		$cmd = "perl $cmdDir/selectRow.pl -f 3 -p -pt \"\" -s  $bed_vs_annot_idpairUniqFile $inBedFile | awk '{print \$2\"\\t\"\$3}' > $outFile";
		print $msgio $cmd, "\n" if $verbose;
		$ret = system ($cmd);
		Carp::croak "$cmd failed: $?\n" if $ret != 0;

		system ("rm -rf  $bed_vs_annot_BedFile $bed_vs_annot_idpairFile $bed_vs_annot_idpairUniqFile") unless $keepCache;
		
	}
	elsif ($summaryMethod eq 'max_num' || $summaryMethod eq 'min_num')
	{
		my $bed_vs_annot_BedFile = "$cache/bed_vs_$name.bed";

		my $cmd = "perl $cmdDir/tagoverlap.pl -d \"/##/\" $bigFlag $verboseFlag $ssFlag $keepCacheFlag --keep-score -region $inBedFile $annotBedFile $bed_vs_annot_BedFile";
		print $msgio $cmd,"\n" if $verbose;

		my $ret = system ($cmd);
		Carp::croak "$cmd failed: $?\n" if $ret != 0;
		my $bed_vs_annot_idpairFile = "$cache/bed_vs_$name.idpair";
																										#interval_id annot_id annot_score
		$cmd = "cut -f 4,5 $bed_vs_annot_BedFile | awk -F \"/##/\" '{print \$1\"\\t\"\$2}' | awk '{print \$2\"\\t\"\$1\"\\t\"\$3}' > $bed_vs_annot_idpairFile";
		print $msgio $cmd, "\n" if $verbose;
		system ($cmd);

		my $bed_vs_annot_idpairUniqFile = "$cache/bed_vs_$name.idpair.uniq";
		$cmd = "perl $cmdDir/uniqRow.pl $verboseFlag -value 2 -c $summaryMethod $bed_vs_annot_idpairFile $bed_vs_annot_idpairUniqFile";
		print $msgio $cmd, "\n" if $verbose;
		$ret = system ($cmd);
		Carp::croak "$cmd failed: $?\n" if $ret != 0;

		#print annotation_id score
		$cmd = "perl $cmdDir/selectRow.pl -f 3 -p -pt \"\" -s  $bed_vs_annot_idpairUniqFile $inBedFile | awk '{print \$2\"\\t\"\$3}' > $outFile";
		print $msgio $cmd, "\n" if $verbose;
		$ret = system ($cmd);
		Carp::croak "$cmd failed: $?\n" if $ret != 0;

		system ("rm -rf  $bed_vs_annot_BedFile $bed_vs_annot_idpairFile $bed_vs_annot_idpairUniqFile");
	
	}
	elsif ($summaryMethod eq 'all')
	{
		my $bed_vs_annot_BedFile = "$cache/bed_vs_$name.bed";

		my $cmd = "perl $cmdDir/tagoverlap.pl -d \"/##/\" $bigFlag $verboseFlag $ssFlag $keepCacheFlag -region $annotBedFile $inBedFile $bed_vs_annot_BedFile";
		print $msgio $cmd,"\n" if $verbose;

		my $ret = system ($cmd);
		Carp::croak "$cmd failed: $?\n" if $ret != 0;
		my $bed_vs_annot_idpairFile = "$cache/bed_vs_$name.idpair";
	
		$cmd = "cut -f 4 $bed_vs_annot_BedFile | awk -F \"/##/\" '{print \$1\"\\t\"\$2}' > $bed_vs_annot_idpairFile";
		print $msgio $cmd, "\n" if $verbose;
		system ($cmd);

		#print all names joined together
		$cmd = "perl $cmdDir/selectRow.pl -f 3 -p -pt \"\" -s  $bed_vs_annot_idpairFile $inBedFile | awk '{print \$2}' > $outFile";
		print $msgio $cmd, "\n" if $verbose;
		$ret = system ($cmd);
		Carp::croak "$cmd failed: $?\n" if $ret != 0;

		system ("rm -rf  $bed_vs_annot_BedFile $bed_vs_annot_idpairFile");
	}
	else
	{
		Carp::croak "unknown summary method: $summaryMethod\n";
	}
}	








