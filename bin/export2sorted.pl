#!/usr/bin/perl -w
# Convert ELAND export.txt to sorted.txt

use strict;
use Time::HiRes qw(time);

if(@ARGV < 2){
	print "Usage: $0 in_export.txt out_sorted.txt\n";
	exit 0;
}
if($ARGV[0] eq '-') {*EXPORT = *STDIN;}
else {open EXPORT, '<', $ARGV[0] or die "Open input export.txt error: $!\n";}
my $unsorted = '.' . time . '.unsorted.tmp';
open UNSORTED, '>', $unsorted or die "Try to open file for output error: $!\n";


# Filter export.txt and put the results in a temp file.
while(<EXPORT>){
	chomp;
	my @cells =	split /\t/;
	my($machine, $run, $lane, $tile, $xcor, $ycor, $index,
		$rn, $read, $qual, 
		$chrom, $contig, $pos, $strand, $desc, 
		$score, $psco, $pchr, $pcont, $poff, $pstnd,
		$filtering) = @cells;
	next if $filtering eq 'N';
	# Remove 'No match', 'QC failed', 'Repeat masked' and 'Multiple hits' reads.
	# Only uniquely aligned reads are retained.
	next if ($chrom eq 'NM' or $chrom eq 'QC' or $chrom eq 'RM'
		or $chrom =~ /^([0-9]+):([0-9]+):([0-9]+)$/);
	pop @cells;	# remove the filtering symbol.
	print UNSORTED join("\t", @cells), "\n";
}
close EXPORT;
close UNSORTED;

# Sort the filtered temp file.
if($ARGV[1] eq '-') {
	`sort -k11,11 -k13,13n -t\"\t\" $unsorted`;
}
else{
	my $sorted = $ARGV[1];
	`sort -k11,11 -k13,13n -t\"\t\" -o $sorted $unsorted`;
}
`rm -f $unsorted`



# An example sorted.txt line split by TAB.
#HWUSI-EAS1713
#0018
#6
#118
#11994
#9213
#0
#1
#ATCTTCTCCTAAGTATCATCCTGAAGAACAAAATTC
#ffefffffffffffffffffffffefeefffefdff
#mm_ref_chr10.fa
#
#3000377
#F
#36
#42
#
#
#
#
#
