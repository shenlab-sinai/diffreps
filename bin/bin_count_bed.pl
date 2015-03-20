#!/usr/bin/perl -w
#
# Bin the genome and count sequence reads for multiple BED files.
# by Li Shen, June 2010.
#


use strict;
use Getopt::Long;
use MyShortRead::MyShortRead;
use MyShortRead::SRBed;

my @bedfiles;	# array of BED files.
my $chrlen;	# chromosome length file.
my $binsize = 1000;	# default bin size = 1kb.
my $fragsize = 0;	# default fragment size = 0, means no read position shift.
if(@ARGV < 1 or !GetOptions('bedfiles=s{1,}' => \@bedfiles,
							'chrlen=s' => \$chrlen,
							'binsize:i' => \$binsize,
							'fragsize:i' => \$fragsize)){
	print "Usage: $0 --bedfiles BED_files --chrlen chromosome_length_txt [--bin_size=1000] [--fragsize=0]\n";
	print "--fragsize	Default=0, means no sequence tag position shift.\n";
	exit 0;
}


# Step0: read chromosome lengths. 
my @chrlen_ordered = read_chrlen_ordered($chrlen);	# ordered array.
my %chrlen_tbl = read_chrlen_tbl($chrlen);	# hash table.

# Step1: create BED file SRBed objects; Separate BED file into trunks of chromosomes.
my @arr_srbed;
foreach my $bed(@bedfiles){
	my $srbed = new MyShortRead::SRBed($bed,$binsize,$fragsize);
	$srbed->create_chr_list(\%chrlen_tbl);
	push @arr_srbed, $srbed;
}

# Step2: for each chromosome, extract bin counts and print.
# Print header line.
print join("\t", "Chrom","Start","End",@bedfiles);
print "\n";
foreach my $r(@chrlen_ordered){
	my $chrom = $r->{chrom};
	# Create bin count vector for current chromosome for all BED files.
	foreach my $b(@arr_srbed) {$b->chr_bin_count($chrom);}
	# Do statistics for each window. All treatments and controls should have the same window_num.
	my $bin_num = $arr_srbed[0]->{chr_list}{$chrom}->{bin_num};
	for my $i(0..$bin_num-1){
		my $start = $i * $binsize + 1;
		my $end = $start + $binsize - 1;
		my @pool_count;
		foreach my $b(@arr_srbed) {push @pool_count, $b->get_bin_count($chrom, $i);}
		print join("\t", $chrom,$start,$end,@pool_count);
		print "\n";
	}
	# Delete bin count vector to save memory. 
	foreach my $r(@arr_srbed) {$r->del_bin_count($chrom);}
}

# Step3: Clean up - delete all temporary chromosome files.
foreach my $r(@arr_srbed) {$r->del_all_chrfiles();}

