#!/usr/bin/perl

use strict;
use warnings;
use Cwd;

##### ##### ##### ##### #####
# Hardwired parameters.

#my $taglen = 4;

##### ##### ##### ##### #####

use Getopt::Std;
use vars qw( $opt_a $opt_c $opt_d $opt_s $opt_p $opt_m $opt_v);

# Usage
my $usage = "
bcsort_fastq_se - Sort barcodes from single-ended fastq files.
		      by
		Brian J. Knaus
		 March 2010

Copyright (c) 2010 Brian J. Knaus.
License is hereby granted for personal, academic, and non-profit use.
Commercial users should contact the author (http://brianknaus.com).

Great effort has been taken to make this software perform its said
task however, this software comes with ABSOLUTELY NO WARRANTY,
not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Usage: perl bcsort_fastq_se.pl options
 required:
  -a	fastq input file.
 optional:
  -c	barcode file name, no quotes, one barcode per line
	[optional, default = barcodes.txt].
  -d	postfix for output directory [optional, default is null
        resulting in outdir = bcsort_qseq_out].
  -s	output fasta file [default is T].
  -m	output fastq file [default is T].
  -v	verbose mode [optional T/F, default is F].

version 4.

";

# command line processing.
getopts('a:c:d:s:m:v:');
die $usage unless ($opt_a);

my ($infa, $bcfile, $outdir, $qseq, $seq, $fq, $verb);

$infa	= $opt_a if $opt_a;
$bcfile	= $opt_c ? $opt_c : "barcodes.txt"; 
$outdir	= $opt_d ? "bcsort_out".$opt_d :"bcsort_out";
$seq	= $opt_s ? $opt_s : "T";
$fq	= $opt_m ? $opt_m : "T";
$verb	= $opt_v ? $opt_v : "F";

##### ##### ##### ##### #####
# Subfunctions,

sub idbc {
	# format id, move barcode to end.
	my @id = split(/\#|\//, shift);
	my $bca = shift;
#	my $fc = shift;
	$id[1] = $bca;
#	return join("", $fc, "_", $id[0], "\#", $id[1], "/", $id[2]);
	my $id = "$id[0]\#$id[1]/$id[2]";
#	return join($id[0], "\#", $id[1], "/", $id[2]);
	return($id);
}

##### ##### ##### ##### #####
# Initialize logfile.

my $bdir = cwd;
mkdir($outdir);
chdir($outdir);

open(my $log, ">", "00_bcsort_log.txt") or die "Can't open my.log: $!";
print $log $usage;
print $log "\n";
print $log "      ----- ----- --*-- ----- -----\n\n";

if ($verb eq "T"){
	print $usage;
  }

# Initial timestmp.
my @stime = localtime(time);
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =
                                                localtime(time);
print $log "Process begun:\t\t";
print $log $mon+1;
print $log "-$mday-";
print $log $year+1900;
print $log " $hour:$min:$sec\n";
print $log "\n      ----- ----- --*-- ----- -----\n\n";

##### ##### ##### ##### #####
# Input barcodes.

chdir($bdir);

# Barcodes
if ($verb eq "T")
  {print "Reading barcodes.\n";}
my $in;

open($in,  "<",  $bcfile)  or die "Can't open $bcfile: $!";

my @temp;
print $log "Barcodes.\n";
my %barcodes;
while (<$in>){
	chomp;
	@temp = split("\t", $_);
	print $log join("\t", @temp, "\n");
	$barcodes{uc($temp[0])} = 0;
}

my $taglen = length($temp[0]);

close $in or die "$in: $!";

print $log "\n      ----- ----- --*-- ----- -----\n\n";


##### ##### ##### ##### #####
# Input data and sort
# into a hash of arrays.
#
# Keys are barcodes.
# Values are arrays of samples.

print $log "Reading in data:\t\t";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =
                                                localtime(time);
print $log $mon+1;
print $log "-$mday-";
print $log $year+1900;
print $log " $hour:$min:$sec\n";
print $log "\n";
close $log or die "$log: $!";

my ($ina, $anom);
open($ina,  "<",  $infa)  or die "Can't open $infa: $!";

chdir($outdir);
open($anom,  ">",  "anomalous_reads.txt")  or die "Can't open anomalous_reads: $!";

my ($ida, $bca, $seqa, $prba, $idb, $bcb, $seqb, $prbb);
my $anom_num = 0;
my $bc_num = 0;
my $ntot = 0;
my %samps;

# Read in.
while (<$ina>){
	$ntot = $ntot + 1;
	# Read a.
	chomp($ida	= substr $_, 1);
	chomp($seqa	= <$ina>);
	$bca	= substr $seqa, 0, $taglen;
	$seqa	= substr $seqa, $taglen;
	chomp($prba	= <$ina>); # Just index.
	chomp($prba	= substr <$ina>, $taglen);

	$ida = idbc($ida, $bca);

	if (exists $barcodes{$bca}){
		# Barcode is expected.
		push @{ $samps{$bca}}, join("\t", $ida, $seqa, $prba); # Hash of arrays.
		$bc_num = $bc_num + 1;
		$barcodes{$bca} = $barcodes{$bca} + 1;
	} else { 
		# Barcode is not expected.
		print $anom "\@$ida\n";
		print $anom "$seqa\n";
		print $anom "+$ida\n";
		print $anom "$prba\n";
		#
		$anom_num = $anom_num + 1;
	}
	if ($ntot % 1000000 == 0){
		open($log, ">>", "00_bcsort_log.txt") or die "Can't open my.log: $!";
		print $log $ntot/1000000, "\tmillion barcodes sorted.\n";
		close $log or die "$log: $!";
	}
	
}

close $ina or die "$ina: $!";
close $anom or die "$anom: $!";

# Report to log.
chdir($outdir);
open($log, ">>", "00_bcsort_log.txt") or die "Can't open my.log: $!";
print $log "\nData read and sorted:\t\t";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =
                                                localtime(time);
print $log $mon+1;
print $log "-$mday-";
print $log $year+1900;
print $log " $hour:$min:$sec\n";

##### ##### ##### ##### #####
# Write data to file.

print $log "\n      ----- ----- --*-- ----- -----\n\n";
print $log "Writing data to file.\n\n";
close $log or die "$log: $!";
my $out1;

# Seq files.
if ($seq eq "T"){
	foreach $bca (keys %samps) {
		open($out1, ">",  $bca."_seq.fa") or die "Can't open output.txt: $!";
		foreach ( @{$samps{$bca}}){
			@temp = split("\t",$_);
			print $out1 ">$temp[0]\n";
			print $out1 "$temp[1]\n";
		}
	close $out1 or die "$out1: $!";
	}
  # Report to log.
  open($log, ">>", "00_bcsort_log.txt") or die "Can't open my.log: $!";
  print $log "seq\tfiles written.\t";
  ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =
	localtime(time);
  print $log $mon+1;
  print $log "-$mday-";
  print $log $year+1900;
  print $log " $hour:$min:$sec\n";
  close $log or die "$log: $!";
}

# Fastq files.
if ($fq eq "T"){
	foreach $bca (keys %samps) {
		open($out1, ">",  $bca."_seq1.fq") or die "Can't open output.txt: $!";
		foreach ( @{$samps{$bca}}){
			@temp = split("\t",$_);
			print $out1 "\@$temp[0]\n";
			print $out1 $temp[1], "\n";
			print $out1 "+$temp[0]\n";
			print $out1 "$temp[2]\n";
		}
		close $out1 or die "out1: $!";
	}
  # Report to log.
  open($log, ">>", "00_bcsort_log.txt") or die "Can't open my.log: $!";
  print $log "fastq\tfiles written.\t";
  ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =
	localtime(time);
  print $log $mon+1;
  print $log "-$mday-";
  print $log $year+1900;
  print $log " $hour:$min:$sec\n";
  close $log or die "$log: $!";
}

open($log, ">>", "00_bcsort_log.txt") or die "Can't open my.log: $!";
print $log "\nData written to file.\t\t";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =
                                                localtime(time);
print $log $mon+1;
print $log "-$mday-";
print $log $year+1900;
print $log " $hour:$min:$sec\n";

##### ##### ##### ##### #####
# Summary.

print $log "\n      ----- ----- --*-- ----- -----\n\n";
print $log "\nTotal reads:\t\t";
print $log $ntot / 1000000, "\tmillion";
print $log "\nBarcoded reads:\t\t";
print $log $bc_num / 1000000, "\tmillion";
print $log "\nAnomalous reads:\t";
print $log $anom_num / 1000000, "\tmillion";


print $log "\n\n";

print $log "Barcode Summary:\n";
while (my($key, $value) = each %barcodes){
	print $log "$key \t=>\t", $value/1000000, "\tmillion", "\n";
}
print "\n";

##### ##### ##### ##### #####
# Finish.

# Reprint start time.
print $log "\n      ----- ----- --*-- ----- -----\n\n";
print $log "Process begun:    ";
print $log $stime[4]+1;
print $log "-$stime[3]-";
print $log $year+1900;
print $log " $stime[2]:$stime[1]:$stime[0]\n";

($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =
                                                localtime(time);

# Log finish time.
print $log "Process finished: ";
print $log $mon+1;
print $log "-$mday-";
print $log $year+1900;
print $log " $hour:$min:$sec\n";
print $log "\n      ----- ----- --*-- ----- -----\n\n";

print $log "bcsort process complete.\n\n";
close $log or die "$log: $!";

##### ##### ##### ##### #####
# EOF.
