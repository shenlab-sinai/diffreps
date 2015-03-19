#!/usr/bin/perl -w
#
# The configuration here only considers chromosome names chrN, chrX, chrY and chrM
# where N is a number.
#

use strict;

my $phred = 0;
if(@ARGV < 1){
  print "Usage: $0 eland_sorted_file [phred_score_cutoff(Default=0)]\n";
  exit 0;
}

if($ARGV[0] eq '-') {*INFILE = *STDIN;}
else {open INFILE, $ARGV[0] or die "Open eland_sorted.txt file error: $!\n";}
$phred = $ARGV[1] if defined $ARGV[1];

my $input = $ARGV[0] eq '-'? 'STDIN' : $ARGV[0];
my $line = 0;
while(<INFILE>){
  chomp;
  my @cells = split /\t/;
  $line++;
  if(@cells < 16){
  	print STDERR "Insufficient info at $input:$line\t$_\n";
	next;
  }
  next if $cells[15] < $phred;	# Skip reads below Phred score cutoff.
  # $cells[10] =~ s/^chr//;
  $cells[10] =~ s/\.fa$//;
  # my $chrom = 'chr' . $cells[10];	# chromosome name.
  my $chrom = $cells[10];	# chromosome name.
  my $rlen = length($cells[8]);	# read length.
  my $strand = $cells[13] eq "F"? "+" : "-";
  if($cells[12] =~ /^[0-9]+$/){	# ensure start position is properly defined.
  	  if($cells[12] > 0){
		  # ELAND: 1-based; BED: 0-based.
		  print join("\t", ($chrom, $cells[12]-1, $cells[12]-1 + $rlen, "N", "0", $strand)), "\n";
	  }else{
	  	  print STDERR "Start position out of range at $input:$line\t$_\n";
		  next;
	  }
  }else{
  	print STDERR "Unproperly defined start position at $input:$line\t$_\n";
	next;
  }
}

close INFILE;

# Example Solexa read.
# HWUSI-EAS1713
# 32
# 1
# 42
# 16389
# 20671
# 0
# 1
# AAGGCCCACCTTGTCTTTACCTTATTTGTTCTAAAT
# dffffffffccfff_dff[cN^_U_fdYfccfadff
# mm_ref_chr10.fa
# 
# 3000111
# F
# 36
# 31
