#!/usr/bin/perl

use strict;
use warnings;
use Cwd;
use POSIX qw(strftime);


##### ##### ##### ##### #####
# Hardwired parameters.

#my $taglen = 4;

##### ##### ##### ##### #####

#Results should be compatible with NCBI submissions (the header modifications currently in the script make the files incompatible with SRA submissions).
#Parameterize the index input to accept either comma separated indexes at the command line or a file of indexes
#Accept multiple indexes of varying length.
#Add a parameter to switch between single and paired end processing (i.e. we would only need one script rather than two).
#Add the ability to recognize different header types, if necessary (different Illumina machines produce different header patterns, which occasionally confuse bcsort).
#- Add the ability to skip blank lines in the fasta file.
#fq or fa

use Getopt::Std;
use vars qw( $opt_a $opt_b $opt_c $opt_d $opt_m $opt_t $opt_v);

# Usage
my $usage = "
bcsort_fastq_pe - Sort barcodes from a pair of fastq files.
		      by
		Brian J. Knaus
		 August 2009

Copyright (c) 2009 Brian J. Knaus.
License is hereby granted for personal, academic, and non-profit use.
Commercial users should contact the author (http://brianknaus.com).

Great effort has been taken to make this software perform its said
task however, this software comes with ABSOLUTELY NO WARRANTY,
not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Usage: perl bcsort_fastq_pe.pl options
 required:
  -a	fasta/q input file a.
  -c	barcodes, either as a comma separated list or a file, one barcode per line, no quotes
  -t  sequence type.  Enter 'fastaa' or 'fastq'
  
 optional:
  -b	fasta/q input file b.
	-m  barcode match type. Default = 4 for PE.
  -d	postfix for output directory. Default = null.
  -v	verbose mode [optional T/F, default is F].

The -m flag specifies what combination of matches the barcodes should have against the sequences.
In the case of single ended reads, there is only one option (a match aginst \"read1\"). The
pattern options for paired end reads are shown below.  If a read matches a barcode it is
shown in caps (i.e. R1 or R2), while no match to a barcode is shown as lowercase (r1 or r2).
1 - R1 AND r2
2 - r1 AND R2
3 - (R1 AND r2) OR (r1 AND R2)
4 - (R1 AND r2) OR (r1 AND R2) OR (R1 AND R2 AND R1=R2) 
5 - R1
6 - R2
7 - R1 OR R2 [R1 value is used if both sequences have a match]
Option 1 is a match only to R1.  Option 2 is a match only to R2.  Option 3 is a match
to either R1 or R2 but not both.  Option 4 is a match to either R1 or R2 or both as long as
they are the identical barcodes.  Option 5 will accept any match to R1, regardless of what, if
anything, matches R2.  Option 6 will accept any match to R2, regardless of what, if
anything, matches R1.  Option 7 will pass any match combination, using the match from R1 if
both reads have matches.

Option 4 accepts the greatest number of reads into non-overlaping partitions, which is why
it is the default.  Using higher option numbers could accept more but if, for instances, the read
files are processed for each barcode separately, the total number of reads across all
partitioned files could be greater than the starting read count.

";

# command line processing.
getopts('a:c:b:f:d:t:v');
die $usage unless ($opt_a);
die $usage unless ($opt_c);

my ($infa, $infb, $fcell, $bcInput, $outdir, $qseq, $seq, $prb, $fq, $verb, $barcodeMatchType);

$infa	= $opt_a;
$infb	= $opt_b;
$bcInput	= $opt_c; 
$outdir	= $opt_d ? "bcsort_out".$opt_d :"bcsort_out";
$verb	= $opt_v ? $opt_v : "F";
if (!$opt_m ) {
	$barcodeMatchType = $opt_m;
}
else {
	$barcodeMatchType = $infb eq '' ? 1 : 4;
}


##### ##### ##### ##### #####
# Initialize logfile.

my $bdir = cwd;
mkdir($outdir);
chdir($outdir);

open(my $log, ">", "00_bcsort_log.txt") or die "Can't open log: $!";
print $log $usage;
print $log "\n";
print $log "      ----- ----- --*-- ----- -----\n\n";

if ($verb eq "T"){
	print $usage;
  }

# Initial timestmp.
my $stime = makeTimestamp();
print $log "Process begun:\t\t".makeTimestamp()."\n";
print $log "\n      ----- ----- --*-- ----- -----\n\n";

##### ##### ##### ##### #####
# Input barcodes.

chdir($bdir);

my $barcodes = getBarcodes($bcInput,$verb,$log);
#my $taglen = length($temp[0]);
print $log "Barcodes used:\n";
print $log join("\n",keys %$barcodes)."\n\n";


print $log "\n      ----- ----- --*-- ----- -----\n\n";


##### ##### ##### ##### #####
# Input data and sort
# into a hash of arrays.
#
# Keys are barcodes.
# Values are arrays of samples.

print $log "Reading in data:\t\t".makeTimestamp()."\n";;

my ($ina, $inb, $anomA, $anomB);
open($ina,  "<",  $infa)  or die "Can't open $infa: $!";
if ($infb){
	open($inb,  "<",  $infb)  or die "Can't open $infb: $!";
}

chdir($outdir);
if ($infb) {
	if ($opt_t eq 'fastq') {
		open($anomA,  ">",  "anomalous_reads_R1.fq")  or die "Can't open anomalous read file: $!";
		open($anomB,  ">",  "anomalous_reads_R2.fq")  or die "Can't open anomalous read file: $!";
	}
	else {
		open($anomA,  ">",  "anomalous_reads_R1.fa")  or die "Can't open anomalous read file: $!";
		open($anomB,  ">",  "anomalous_reads_R2.fa")  or die "Can't open anomalous read file: $!";
	}
}
else {
	if ($opt_t eq 'fastq') {
		open($anomA,  ">",  "anomalous_reads.fq")  or die "Can't open anomalous read file: $!";
	}
	else {
		open($anomA,  ">",  "anomalous_reads.fa")  or die "Can't open anomalous read file: $!";
	}
}

my ($ida, $bca, $seqa, $prba, $idb, $bcb, $seqb, $prbb);
my $anom_num = 0;
my $bc_num = 0;
my $ntot = 0;
my %samps;



# Loop through sequence file(s).
while (<$ina>){
	$ntot = $ntot + 1;
	
	chomp($ida	= $_);
	chomp($seqa	= <$ina>);
	
	if ($infb) {
		chomp($idb	= substr <$inb>, 1);
		chomp($seqb	= <$inb>);
	}
	
	if ($opt_t eq 'fastq') {
		# Read a.
		$prba	= <$ina>; # skip the header
		chomp($prba	= <$ina>); # qualstring.

	# Read b.
		if ($infb) {
			$prbb = <$inb>; #skip the header
			chomp($prbb	= <$inb>); # qualstring
		}
	}
	
	#check for barcode match
	my $filterValA = matchBarcodes($ida,$seqa,$barcodes);
	my $filterValB = "";
	if ($infb) {
		$filterValB = matchBarcodes($idb,$seqb,$barcodes);
	}
	print "A=$filterValA, B=$filterValB\n";
	
	my $printFlag = 0;
	if ($infb) {
		if (($filterValA ne '' and $filterValB eq '') or
				($filterValA eq '' and $filterValB ne '') or
				($filterValA eq $filterValB and $filterValA ne '')) {
			$printFlag = 1;
		}
	}
	else {
		$printFlag = 1 if $filterValA ne '';
	}
	
	#if filters not met, print to anomalous files and go to next sequence
	if ($printFlag == 0){
		print $anomA $ida."\n";
		print $anomA $seqa."\n";
		if ($opt_t eq 'fastq') {
			print $anomA "+\n$prba\n";
		}
		if ($infb){
			print $anomB $idb."\n";
			print $anomB $seqb."\n";
			if ($opt_t eq 'fastq') {
				print $anomB "+\n$prbb\n";
			}
		}
		
		
		
		next;
	}
	
	my $barcode = $filterValA eq '' ? $filterValB : $filterValA;
	print "barcode = $barcode, printflag = $printFlag\n";
	
	#open file handlers for output
	my ($outA,$outB);
	if ($opt_t eq 'fastq') {
		if ($infb ne '') {
			open ($outA,">>",$barcode."_R1.fq");
			open ($outB,">>",$barcode."_R2.fq");
		}
		else {
			open ($outA,">>",$barcode."_seq.fq");
		}
	}
	else {
		if ($infb ne '') {
			open ($outA,">>",$barcode."_R1.fa");
			open ($outB,">>",$barcode."_R2.fa");
		}
		else {
			open ($outA,">>",$barcode."_seq.fa");	
		}
		
		
	}
	
	#output sequences
	print $outA "$ida\n";
	print $outA substr($seqa,length($filterValA))."\n";
	if ($opt_t eq 'fastq') {
		print $outA "+\n";
		print $outA substr($prba,length($filterValA))."\n";
	}
	if ($infb ne '') {
		print $outB "$idb\n";
		print $outB substr($seqb,length($filterValA))."\n";
		if ($opt_t eq 'fastq') {
			print $outB "+\n";
			print $outB substr($prbb,length($filterValA))."\n";
		}
	}

	close $outA;
	if ($infb) {
		close $outB;
	}
	
	
}

close $ina or die "$ina: $!";
if ($infb) {
	close $inb or die "$inb: $!";	
}

#close $anom or die "$anom: $!";
close $anomA;
if ($infb) {
	close $anomB;
}

# Report to log.
chdir($outdir);
#open($log, ">>", "00_bcsort_log.txt") or die "Can't open my.log: $!";
print $log "\nData read and sorted:\t\t".makeTimestamp()."\n";

#chdir($outdir);


=begin
##### ##### ##### ##### #####
# Write data to file.

print $log "\n      ----- ----- --*-- ----- -----\n\n";
print $log "Writing data to file.\n\n";
close $log or die "$log: $!";
my $out1;
#my @temp;

# Seq files.
if ($seq eq "T"){
	foreach $bca (keys %samps) {
		open($out1, ">",  $bca."_seq.fa") or die "Can't open output.txt: $!";
		foreach ( @{$samps{$bca}}){
			@temp = split("\t",$_);
			print $out1 ">$temp[0]\n";
			print $out1 "$temp[1]\n";
			print $out1 ">$temp[3]\n";
			print $out1 "$temp[4]\n";
		}
	close $out1 or die "$out1: $!";
	}
}
# Report to log.
open($log, ">>", "00_bcsort_log.txt") or die "Can't open my.log: $!";
print $log "seq\tfiles written (if applicable).\t";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =
	localtime(time);
print $log $mon+1;
print $log "-$mday-";
print $log $year+1900;
print $log " $hour:$min:$sec\n";
close $log or die "$log: $!";

# Prb files.
if ($prb eq "T"){
	foreach $bca (keys %samps) {
		open($out1, ">",  $bca."_prb.fa") or die "Can't open output.txt: $!";
		foreach ( @{$samps{$bca}}){
			@temp = split("\t",$_);
			print $out1 ">$temp[0]\n";
			print $out1 $temp[2], "\n";
#			print $out1 asc2phred2($temp[2]), "\n";
			print $out1 ">$temp[3]\n";
			print $out1 $temp[5], "\n";
#			print $out1 asc2phred2($temp[5]), "\n";
		}
		close $out1 or die "$out1: $!";
	}
}
# Report to log.
open($log, ">>", "00_bcsort_log.txt") or die "Can't open my.log: $!";
print $log "prb\tfiles written (if applicable).\t";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =
	localtime(time);
print $log $mon+1;
print $log "-$mday-";
print $log $year+1900;
print $log " $hour:$min:$sec\n";
close $log or die "$log: $!";

# Fastq files.
if ($fq eq "T"){
	foreach $bca (keys %samps) {
		open($out1, ">",  $bca."_seq1.fq") or die "Can't open output.txt: $!";
		open(my $out2, ">",  $bca."_seq2.fq") or die "Can't open output.txt: $!";
		foreach ( @{$samps{$bca}}){
			@temp = split("\t",$_);
			print $out1 "\@$temp[0]\n";
			print $out1 $temp[1], "\n";
			print $out1 "+$temp[0]\n";
			print $out1 "$temp[2]\n";
			print $out2 "\@$temp[3]\n";
			print $out2 $temp[4], "\n";
			print $out2 "+$temp[3]\n";
			print $out2 $temp[5], "\n";
		}
		close $out1 or die "out1: $!";
		close $out2 or die "out2: $!";
	}
}
# Report to log.
open($log, ">>", "00_bcsort_log.txt") or die "Can't open my.log: $!";
print $log "fastq\tfiles written (if applicable).\t";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =
	localtime(time);
print $log $mon+1;
print $log "-$mday-";
print $log $year+1900;
print $log " $hour:$min:$sec\n";
close $log or die "$log: $!";

open($log, ">>", "00_bcsort_log.txt") or die "Can't open my.log: $!";
print $log "\nData written to file.\t\t";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =
                                                localtime(time);
print $log $mon+1;
print $log "-$mday-";
print $log $year+1900;
print $log " $hour:$min:$sec\n";
=cut
##### ##### ##### ##### #####
# Summary.

print $log "\n      ----- ----- --*-- ----- -----\n\n";
print $log "\nTotal reads:\t\t";
print $log $ntot / 1000000, "\tmillion";
print $log "\nBarcoded reads:\t\t";
print $log $bc_num / 1000000, "\tmillion";
print $log "\nNon-barcoded reads:\t";
print $log ($ntot - $bc_num) / 1000000, "\tmillion";
print $log "\nAnomalous reads:\t";
print $log $anom_num / 1000000, "\tmillion";


print $log "\n\n";

print $log "Barcode Summary:\n";
while (my($key, $value) = each %$barcodes){
	print $log "$key \t=>\t", $value/1000000, "\tmillion", "\n";
}
print "\n";

##### ##### ##### ##### #####
# Finish.

# Reprint start time.
print $log "\n      ----- ----- --*-- ----- -----\n\n";
print $log "Process begun: $stime ";

# Log finish time.
print $log "Process finished: ".makeTimestamp()."\n\n";
print $log "\n      ----- ----- --*-- ----- -----\n\n";

print $log "bcsort process complete.\n\n";
close $log or die "$log: $!";


##### ##### ##### ##### #####
# Subfunctions,

sub asc2phred {
  my $asc = shift;
  my $phred = ord(substr($asc, 0, 1)) - 64;
  my $index = 1;
  while ($index < length($asc))
    {
    $phred = join " ", $phred, ord(substr($asc, $index, 1)) - 64;
    $index = $index +1;
    }
  return $phred;
  }

sub asc2phred2 {
	my @temp = split(//, $_);
	foreach(@temp){
		$_ = ord($_) - 64;
	}
	return(join " ", @temp);
}

#sub qseq2id {
#	my @temp = split("\t", $_);
#	return join(":", $temp[0], $temp[2], $temp[3], $temp[4], $temp[5] );
#}




sub matchBarcodes {
	my $header = shift;
	my $seq = shift;
	my $barcodes = shift;
	my $matchCount = 0;
	my $matchBarcode = "";
	foreach my $bc (keys %$barcodes) {
		#print "bc=$bc, seq = ".substr($seq,0,10)."\n";
		if ($seq =~ /^$bc/) {
			$matchBarcode = $bc;
			$matchCount++;
		}
	}
	if ($matchCount>1) {
		print STDERR "$header has matched more than one barcode.  This is very very bad, I'm shutting down for an oil bath."
	}
	return $matchBarcode;
}


sub idbc {
	# id tab bc as input.
	my @id = split(/\#|\//, shift);
	my $bca = shift;
	my $fc = shift;
	$id[1] = $bca;
	my $returnstring = join("", $fc, "_", $id[0], "\#", $id[1]);
	if ($id[2]) {
		$returnstring .= "/$id[2]"
	}
	return $returnstring;
}

sub getBarcodes {
	my $bcInput = shift;
	my $verb = shift;
	my $logFH = shift;
	
	
	if ($verb eq "T")
		{print "Reading barcodes.\n";}
	
	my %barcodes;
	if (-f $bcInput) {
		
		my $in;
		open($in,  "<",  $bcInput)  or die "Can't open $bcInput: $!";
		
		my @temp;
		print $logFH "Barcodes.\n";
		while (<$in>){
			chomp;
			@temp = split("\t", $_);
			$barcodes{uc($temp[0])} = 0;
		}
	}
	
	else {
		foreach my $bc (split(",",$bcInput)) {
			$bc =~ s/^\s+|\s+$//g;
			$barcodes{uc($bc)} = 0;
		}
	}
	
	print $logFH "Barcodes.\n";
	print $logFH join(",",keys(%barcodes))."\n";
	
	return \%barcodes;
}


sub makeTimestamp {
	return strftime "%d-%b-%Y %I:%M:%S", localtime;
}

##### ##### ##### ##### #####
# EOF.
