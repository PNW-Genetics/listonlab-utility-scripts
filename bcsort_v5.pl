#!/usr/bin/perl

use strict;
use warnings;
use Cwd;
use POSIX qw(strftime);
use IO::File;
#use Time::HiRes qw(gettimeofday tv_interval)

#print "this\n";
#=begin
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
use vars qw( $opt_a $opt_b $opt_c $opt_m $opt_t $opt_o);

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
  -a  fasta/q input file a.
  -c  barcodes, either as a comma separated list or a file, one barcode per line, no quotes
  -t  sequence type.  Enter 'fastaa' or 'fastq'
  
 optional:
  -b  fasta/q input file b.
  -m  barcode match type. Default = 4 for PE.
  -o  Output directory name (will be created if doesnt exist).  Default = bcsort_out.

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
getopts('a:c:b:f:t:m:o:v');
die $usage unless ($opt_a);
die $usage unless ($opt_c);
die $usage unless ($opt_t);

my ($infa, $infb, $fcell, $bcInput, $outdir, $qseq, $seq, $prb, $fq, $barcodeMatchType);

$infa	= $opt_a;
$infb	= $opt_b;
$bcInput	= $opt_c;
$outdir	= $opt_o || "bcsort_out";
if ($opt_m) {
	$barcodeMatchType = int($opt_m);
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

print $log $usage."\n";
print $log "      ----- ----- --*-- ----- -----\n\n";

print $log "Command run:\n";
print $log "$0 -a $opt_a ";
print $log "-b $opt_b " if $opt_b;
print $log "-c $opt_c -t $opt_t -m $barcodeMatchType ";
print $log "-o $outdir ";
print $log "\n\n";
print $log "      ----- ----- --*-- ----- -----\n\n";

# Initial timestmp.
my $stime = makeTimestamp();
print $log "Process begun:\t\t".makeTimestamp()."\n";
print $log "\n      ----- ----- --*-- ----- -----\n\n";

##### ##### ##### ##### #####
# Input barcodes.

chdir($bdir);

my $barcodes = getBarcodes($bcInput,$log);
#my $taglen = length($temp[0]);


print $log "\n      ----- ----- --*-- ----- -----\n\n";



##### ##### ##### ##### #####
# Input data and sort
# into a hash of arrays.
#
# Keys are barcodes.
# Values are arrays of samples.

print $log "Reading in data:\t\t".makeTimestamp()."\n";;

#my ($ina, $inb, $anomA, $anomB);

#my ($ida, $seqa, $prba, $idb, $seqb, $prbb);
my $anom_num = 0;
my $bc_num = 0;
my $ntot = 0;
my %samps;

#process paired end reads
if ($infb) {
	
	#make output files for each barcode:
	chdir($outdir);
	my %fileHandles;
	foreach my $bc (keys %$barcodes) {
		if ($opt_t eq 'fastq'){
			$fileHandles{$bc}->{'R1'}=IO::File->new(">$bc\_R1.fq");
			$fileHandles{$bc}->{'R2'}=IO::File->new(">$bc\_R2.fq");
		}
		else{
			$fileHandles{$bc}->{'R1'}=IO::File->new(">$bc\_R1.fa");
			$fileHandles{$bc}->{'R2'}=IO::File->new(">$bc\_R2.fa");
		}
	}

	
	# Assign an anonymous subroutine to check if the match criteria are met
	# Only need to do this once per run, so better to do it out here than in input loop
	my $passesMatch;
	if ($barcodeMatchType == 1) {
		$passesMatch = sub{
			my $R1 = shift;
			my $R2 = shift;
			if ($R1 ne '' && $R2 eq '') {
				return 1;
			}
			else {
				return 0;
			}
		}
	}
	elsif ($barcodeMatchType == 2) {
		$passesMatch = sub{
			my $R1 = shift;
			my $R2 = shift;
			if ($R1 eq '' && $R2 ne '') {
				return 1;
			}
			else {
				return 0;
			}
		}
	}
	elsif ($barcodeMatchType == 3) {
		$passesMatch = sub{
			my $R1 = shift;
			my $R2 = shift;
			if (($R1 eq '' && $R2 ne '') ||
					($R1 ne '' && $R2 eq '')) {
				return 1;
			}
			else {
				return 0;
			}
		}
	}
	elsif ($barcodeMatchType == 4) {
		$passesMatch = sub{
			my $R1 = shift;
			my $R2 = shift;
			if (($R1 eq '' && $R2 ne '') ||
					($R1 ne '' && $R2 eq '') ||
					($R1 eq $R2 && $R1 ne '')) {
				return 1;
			}
			else {
				return 0;
			}
		}
	}
	elsif ($barcodeMatchType == 5) {
		$passesMatch = sub{
			my $R1 = shift;
			my $R2 = shift;
			if ($R1 ne '') {
				return 1;
			}
			else {
				return 0;
			}
		}
	}
	elsif ($barcodeMatchType == 6) {
		$passesMatch = sub{
			my $R1 = shift;
			my $R2 = shift;
			if ($R2 ne '') {
				return 1;
			}
			else {
				return 0;
			}
		}
	}
	elsif ($barcodeMatchType == 7) {
		$passesMatch = sub{
			my $R1 = shift;
			my $R2 = shift;
			if ($R1 ne '' || $R2 ne '') {
				return 1;
			}
			else {
				return 0;
			}
		}
	}

	chdir($bdir);
	my ($ina, $inb, $anomA, $anomB);
	open($ina,  "<",  $infa)  or die "Can't open $infa: $!";
	open($inb,  "<",  $infb)  or die "Can't open $infb: $!";


	chdir($outdir);	
	if ($opt_t eq 'fastq') {
		open($anomA,  ">",  "anomalous_reads_R1.fq")  or die "Can't open anomalous read file: $!";
		open($anomB,  ">",  "anomalous_reads_R2.fq")  or die "Can't open anomalous read file: $!";
	}
	else {
		open($anomA,  ">",  "anomalous_reads_R1.fa")  or die "Can't open anomalous read file: $!";
		open($anomB,  ">",  "anomalous_reads_R2.fa")  or die "Can't open anomalous read file: $!";
	}

	#set paramaters based on fasta or fastq
	my ($outSuffix);
	
	if ($opt_t eq 'fastq'){
		$outSuffix = "_seq.fq";
	}
	else {
		$outSuffix = "_seq.fa";
	}
	
	my ($ida, $idb, $seqa, $seqb, $prba, $prbb);
	# Loop through sequence file(s).
	while (my $inLine = <$ina>){
		$ntot = $ntot + 1;
		
		

		#read the sequence
		chomp($ida	= $inLine);
		chomp($seqa	= <$ina>);

		chomp($idb	= <$inb>);
		chomp($seqb	= <$inb>);
		
		if ($opt_t eq 'fastq') {
			$prba	= <$ina>; # skip the header
			chomp($prba	= <$ina>); # qualstring
			$prbb = <$inb>; #skip the header
			chomp($prbb	= <$inb>); # qualstring
		}
		
		#check for barcode match

		
		my $filterValA = matchBarcodes($ida,$seqa,$barcodes);
		my $filterValB = matchBarcodes($idb,$seqb,$barcodes);
		
		#print "A=$filterValA, B=$filterValB\n";
		
		my $printFlag = $passesMatch->($filterValA,$filterValB);
		#print "printflag=$printFlag\n";	
		
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
		}
			
		else {
			$bc_num++;
			my $barcode = $filterValA eq '' ? $filterValB : $filterValA;
			$barcodes->{$barcode}+=1;
			#print "barcode = $barcode, printflag = $printFlag\n";
			
			#open file handlers for output
			
			my $outA = $fileHandles{$barcode}->{'R1'};
			my $outB = $fileHandles{$barcode}->{'R2'};
			
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
			
			
		}
		
		#print $ntot."\n";
		if ($ntot % 500000 == 0){
			print "Finished processing $ntot lines: ".makeTimestamp()."\n";
		}
	}

	
	close $ina;
	close $inb;	
	close $anomA;
	close $anomB;


}

#process SE reads
else {

}


# Report to log.
chdir($outdir);
#open($log, ">>", "00_bcsort_log.txt") or die "Can't open my.log: $!";
print $log "\nData read and parsed:\t\t".makeTimestamp()."\n";

#chdir($outdir);


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
print $log "Process begun: $stime\n";

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
	my $logFH = shift;
	
	my %barcodes;
	if (-f $bcInput) {
		
		my $in;
		open($in,  "<",  $bcInput)  or die "Can't open $bcInput: $!";
		
		my @temp;
		while (<$in>){
			next if $_=~/^\n+$/;
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
#=cut
##### ##### ##### ##### #####
# EOF.
