#!/usr/bin/perl -w
#
# Identify sliding windows above read count cutoff and output continuous
# regions into a BED file.
#

use strict;
use Getopt::Long;
use MyShortRead::MyShortRead;
use MyBioinfo::Common;
use MyShortRead::SRBed;
use diffReps::ChromaModSite;

my @treatment;	# array of treatment file names.
my @control;	# array of control file names.
my($nsd,$step,$window,$chrlen,$fragsize,$bed);
$nsd = 2;			# default number of std above mean as cutoff.
$step = 100;		# default step/bin size.
$window = 1000;
$fragsize = 200;	# default fragment size.
if(@ARGV < 1 or !GetOptions('treatment=s{1,}' => \@treatment,
							'control=s{1,}' => \@control,
							'chrlen=s' => \$chrlen,
							'bed=s' => \$bed,
							'nsd:i' => \$nsd,
							'window:i' => \$window,
							'step:i' => \$step,
							'fragsize:i' => \$fragsize)){
	print "\n\$\$ Output Window Coverage to A BED File. \$\$\n\n";
	print "Usage: $0 --treatment bed_file1...[bed_fileN] --control bed_file1...[bed_fileN] --chrlen chrom_length_file --bed output_BED\n\n";
	print "Optional parameters(defaults in parentheses):\n";
	print "--nsd($nsd)                 Number of std above mean as count cutoff.\n";
	print "--window($window)           Window size.\n";
	print "--step($step)              Window moving step size.\n";
	print "--frag($fragsize)              ChIP-seq library fragment size. Use to shift read positions.\n\n";
	exit 0;
}

# Step0: read chromosome length.
my @chrlen_ordered = read_chrlen_ordered($chrlen);	# ordered array.
my %chrlen_tbl = read_chrlen_tbl($chrlen);	# hash table.

# Step1: create treatment and control SRBed objects;
# separate BED file into trunks of chromosomes.
my @tr_srbed;
my @co_srbed;
foreach my $bedfile(@treatment){
	my $srbed = new MyShortRead::SRBed($bedfile,$window,$fragsize);
	$srbed->create_chr_list(\%chrlen_tbl);
	push @tr_srbed, $srbed;
}
foreach my $bedfile(@control){
	my $srbed = new MyShortRead::SRBed($bedfile,$window,$fragsize);
	$srbed->create_chr_list(\%chrlen_tbl);
	push @co_srbed, $srbed;
}

# Step2: find mean, std and normalization constants for all samples.
my @all_beds = (@tr_srbed, @co_srbed);
my $ms_res = mean_std(\@all_beds, \@chrlen_ordered);
my(@tr_mean, @co_mean, @tr_std, @co_std);
split_array($ms_res->{mean}, scalar @tr_srbed, \@tr_mean, \@co_mean);
split_array($ms_res->{std}, scalar @tr_srbed, \@tr_std, \@co_std);
# Reset bin and window sizes for sliding window procedure.
foreach my $b(@tr_srbed){
	$b->set_bin_size($step);
	$b->set_window_size($window);
}
foreach my $b(@co_srbed){
	$b->set_bin_size($step);
	$b->set_window_size($window);
}

# Step3: for each chromosome, do statistics on windows and merge continuous windows into one region.
open BED, '>', $bed or die "Open BED file error: $!\n";
my $global_data = {	# hash to hold global data.
	bedTr => \@tr_srbed,	# array of references to treatment SRBed objects.
	bedCo => \@co_srbed,	# array of references to control SRBed objects.
	meanTr => \@tr_mean,	# mean bin count for treatment.
	meanCo => \@co_mean,	# mean bin count for control.
	stdTr => \@tr_std,		# bin count std for treatment.
	stdCo => \@co_std,		# bin count std for control.
	nsd => $nsd,			# number of std above mean as count cutoff.
	winSize => $window,		# window size.
	step => $step,			# moving step size.
};
my $slideWindow = new diffReps::SlideWindow;	# sliding window.
foreach my $r(@chrlen_ordered){
	my $chrom = $r->{chrom};
	# Create bin count and window count vectors for current chromosome for treatment and control.
	my $window_num;	# should be the same for all treatment and control samples.
	$window_num = create_chr_cntvec(\@tr_srbed, $chrom) or 
		die "Create count vector failed due to zero treatment sample!";
	$window_num = create_chr_cntvec(\@co_srbed, $chrom) or 
		die "Create count vector failed due to zero control sample!";
	# Prepare window to slide.
	$slideWindow->set_chrom($global_data, $chrom, $window_num, 1);
	my $wincovRegion = new WinCovRegion;	# class to store window coverage region.

	# Search for significant windows and then expand the region.
	while($slideWindow->get_winN()>=0){
		if($slideWindow->isAboveCutoff($global_data)){
			if($wincovRegion->isContinuous($global_data, $slideWindow)){
				$wincovRegion->expand($global_data, $slideWindow);
			}
			else{
				$wincovRegion->print_bed(\*BED);
				$wincovRegion->init($global_data, $slideWindow);
			}
		}
		$slideWindow->move_forward($global_data, 1);
	}
	# Print the last significant region left in memory.
	$wincovRegion->print_bed(\*BED);
	
	# Delete bin count and window count vectors to save memory. 
	foreach my $r(@tr_srbed){
		$r->del_bin_count($chrom);
		$r->del_win_count($chrom);
	}
	foreach my $r(@co_srbed){
		$r->del_bin_count($chrom);
		$r->del_win_count($chrom);
	}
}
close BED;

# Step4: Clean up - delete all temporary chromosome files.
foreach my $r(@tr_srbed) {$r->del_all_chrfiles();}
foreach my $r(@co_srbed) {$r->del_all_chrfiles();}


################ START Subroutines ###################
# Create bin and window count vectors for a chromosome.
# Return number of windows.
sub create_chr_cntvec{
	my($rbeds, $chrom) = @_;
	if(!@{$rbeds}) {return 0;}
	foreach my $b(@{$rbeds}){
		$b->chr_bin_count($chrom);
		$b->chr_win_count($chrom);
	}
	return $rbeds->[0]->{chr_list}{$chrom}->{window_num};
}

# Go through a set of SRBed objects and calculate the mean, std and normalization constants.
sub mean_std{
	my($rbeds,$chrlen) = @_;
	my(@sums,@sumsqs);	# array of sum and sum of squares.
	foreach my $b(@{$rbeds}) {	# initialize arrays.
		push @sums, 0;
		push @sumsqs, 0;
	}
	my $res = {	# result hash table.
		mean =>	[],	# array of means.
		std =>	[],	# array of stds.
	};

	#### Calculate mean and std. ####
	my $nbin = 0;	# number of bins in the genome.
	foreach my $r(@{$chrlen}){
		my $chrom = $r->{chrom};
		foreach my $b(@{$rbeds}) {$b->chr_bin_count($chrom);}
		my $maxN = $rbeds->[0]->{chr_list}{$chrom}->{bin_num};	# num of bins in chrom.
		$nbin += $maxN;
		for my $i(0..$maxN-1){	# go through each bin.
			my @cnt;
			my $j = 0;	# iterator for SRBed objects.
			foreach my $b(@{$rbeds}) {
				push @cnt, $b->get_bin_count($chrom, $i);
				if($cnt[$j] > 0){
					$sums[$j] += $cnt[$j];
					$sumsqs[$j] += $cnt[$j]**2;
				}
				$j++;
			}
		}
		# Delete bin count vectors to free memory.
		foreach my $b(@{$rbeds}) {$b->del_bin_count($chrom);}
	}
	foreach my $s(@sums) {push @{$res->{mean}}, $s/$nbin;}
	my $j = 0;	# iterator for SRBed objects.
	foreach my $s(@sumsqs) {push @{$res->{std}}, sqrt($s/$nbin - $res->{mean}[$j++]**2);}

	return $res;
}

# Split an array into two halves.
sub split_array{
	my($O,$len,$A,$B) = @_;
	@{$B} = @{$O};
	@{$A} = splice @{$B}, 0, $len;
}

################## END Subroutines ####################


# A class for determining regions of continuous windows above count cutoff.
# The regions can be output as BED records.
package WinCovRegion;
use strict;

# Constructor.
sub new{
	my($class) = @_;
	my $self = {
		chrom => undef,
		start => undef,
		end => undef
	};
	bless $self, $class;
	return $self;
}

# Initialize a region using a window.
sub init{
	my($self,$gdt,$rwin) = @_;
	$self->{chrom} = $rwin->get_chrom();
	$self->{start} = $rwin->get_winN()*$gdt->{step} + 1;
	$self->{end} = $rwin->get_winN()*$gdt->{step} + $gdt->{winSize};
}

# Is the window continuous (overlap or adjacent) to the existing region?
sub isContinuous{
	my($self,$gdt,$rwin) = @_;
	return 0 if !defined $self->{chrom};
	my $winSta = $rwin->get_winN()*$gdt->{step} + 1;
	return ($winSta - $self->{end} <= 1);
}

# Expand existing region with new coordinates.
sub expand{
	my($self,$gdt,$rwin) = @_;
	my $winEnd = $rwin->get_winN()*$gdt->{step} + $gdt->{winSize};
	$self->{end} = $winEnd;
}

# Return coordinates of the existing region.
sub get_coord{
	my($self) = @_;
	return {
		chrom => $self->{chrom},
		start => $self->{start},
		end => $self->{end}
	};
}

# Print current region coordinates as a BED record.
sub print_bed{
	my($self,$h) = @_;
	if(defined $self->{chrom}){
		print $h join("\t", $self->{chrom}, $self->{start}, $self->{end}), "\n";
	}
}





