#!/usr/bin/perl -w
# Remove redundant reads from sorted BED file.

use strict;

## Program parameters ##
if(@ARGV < 3){
	print "Usage: $0 in_sorted.bed out_nonredundant.bed dup_allowed\n";
	exit 0;
}
if($ARGV[0] eq '-') {*SORTED = *STDIN;}
else {open SORTED, '<', $ARGV[0] or die "Open input sorted.bed error:$!\n";}
if($ARGV[1] eq '-') {*NONR = *STDOUT;}
else {open NONR, '>', $ARGV[1] or die "Open output nonr.bed error:$!\n";}
my $dup_allowed = $ARGV[2];
die "\$dup_allowed must be no less than 1! Exit.\n" if $dup_allowed < 1;

## Iterate through BED file ##
my $cnt_pos = 0;	# counter for positive strand.
my $cnt_neg = 0;	# counter for negative strand.
my $cur_chrom = '';
my $cur_start = -1;
while(<SORTED>){
	chomp;
	my($chrom, $start, $end, $name, $score, $strand) = split /\t/;
	if($chrom ne $cur_chrom or $start != $cur_start){
		$cur_chrom = $chrom;
		$cur_start = $start;
		$cnt_pos = 0;
		$cnt_neg = 0;
	}
	if($strand eq '+'){
		$cnt_pos++;
		next if $cnt_pos > $dup_allowed;
	}
	else {
		$cnt_neg++;
		next if $cnt_neg > $dup_allowed;
	}
	print NONR "$chrom\t$start\t$end\t$name\t$score\t$strand\n";
}
close SORTED;
close NONR;
