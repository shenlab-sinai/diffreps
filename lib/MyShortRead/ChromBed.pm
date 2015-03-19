#
# ChromBed - A class for extracting short read information for one chromosome.
#

use strict;
use 5.006;
package MyShortRead::ChromBed;
use POSIX qw(ceil);
use MyShortRead::MyShortRead qw(bin_genome_count bin_genome_count_direction);
use MyBioinfo::Common qw(mean median);
our $VERSION = '1.00';

# Class constructor.
sub new{
	my $class = shift;
	my $self = {};
	$self->{chrom_name} = undef;
	$self->{chrom_len} = undef;
	$self->{file_name} = undef;	# chromosome file name.
	$self->{bin_num} = undef;	# number of bins for the chromosome.
	$self->{bin_count} = [];	# vector for bin read count.
	$self->{bin_count_F} = [];	# bin Forward strand read count.
	$self->{bin_count_R} = [];	# bin Reverse strand read count.
	$self->{window_num} = undef;	# number of sliding windows.
	$self->{window_count} = [];	# vector for window read count.
	$self->{window_count_F} = [];	# window read count for Forward strand.
	$self->{window_count_R} = [];	# window read count for Reverse strand.
	bless($self, $class);

	# User supplied chromosome name, length and file name.
	if(@_ >= 3){
		$self->{chrom_name} = $_[0];
		$self->{chrom_len} = $_[1];
		$self->{file_name} = $_[2];
	}

	return $self;
}

# Set chromosome information: name, length and file name.
# Clean memory data.
sub set_chr_info{
	my $self = shift;
	if(@_ >= 3){
		$self->{chrom_name} = $_[0];
		$self->{chrom_len} = $_[1];
		$self->{file_name} = $_[2];
		$self->{bin_num} = undef;
		$self->{bin_count} = [];
	}
}

# Create bin count vector for the specified chromosome file.
sub do_bin_count{
	my $self = shift;
	if(@_ < 2) {warn "Bin or fragment size not specified. Bin count not performed. Exit.\n"; return;}
	if(!(defined $self->{chrom_name} and defined $self->{chrom_len} and 
		defined $self->{file_name})){
		warn "Chromosome information missing! Cannot continue. Exit.\n";
		return;
	}
	my($bin_size,$frag_size) = @_;
	# Read bin size and calculate bin number.
	$self->{bin_num} = ceil($self->{chrom_len} / $bin_size);
	# Create bin count vector.
	bin_genome_count($self->{file_name},$bin_size,$self->{chrom_len},
		$frag_size,$self->{bin_count});
}

# Create bin count with direction for the specified chromosome file.
sub do_bin_count_direction{
	my $self = shift;
	if(@_ < 2) {warn "Bin or fragment size not specified. Bin count not performed. Exit.\n"; return;}
	if(!(defined $self->{chrom_name} and defined $self->{chrom_len} and defined $self->{file_name})){
		warn "Chromosome information missing! Cannot continue. Exit.\n";
		return;
	}
	my($bin_size,$frag_size) = @_;
	# Read bin size and calculate bin number.
	$self->{bin_num} = ceil($self->{chrom_len} / $bin_size);
	# Create bin count vector.
	bin_genome_count_direction($self->{file_name},$bin_size,$self->{chrom_len},
		$frag_size,$self->{bin_count_F},$self->{bin_count_R});
}

# Smooth bin count signals by applying window functions.
sub smooth_bin_count{
	my $self = shift;
	if(@_ < 2) {warn "Smoothing width or function not specified. Exit.\n"; return;}
	if(!defined $self->{bin_count}){
		warn "Bin count vector not exists. Do bin count first! Exit.\n";
		return;
	}
	my($smooth_width,$fun_str) = @_;
	if($smooth_width <= 1){
		warn "Smooth width must be larger than one. Exit.\n";
		return;
	}
	if($smooth_width % 2 != 1){
		warn "Smooth width must be odd number. Convert it to nearest odd number.\n";
		$smooth_width++;
	}
	my $half = ($smooth_width-1)/2;
	if($self->{bin_num} <= $half){
		warn "Bin vector too small to calculate smoothed signals.\n";
		return;
	}
	my @mov_arr;	# Array containing bin counts of moving window.
	# Build the initial moving array.
	for my $i(0..$half) {push @mov_arr, $self->{bin_count}[$i];}
	if($fun_str eq 'mean') {$self->{bin_count}[0] = mean(@mov_arr);}
	elsif($fun_str eq 'median') {$self->{bin_count}[0] = median(@mov_arr);}
	# Iteratively updating moving array and calculate smoothed signals.
	for(my $i = 1; $i < $self->{bin_num}; $i++){
		if($i - $half > 0) {shift @mov_arr;}	# Remove the left element.
		if($i + $half < $self->{bin_num}) {push @mov_arr, $self->{bin_count}[$i+$half];}	# Add new element to the right.
		if($fun_str eq 'mean') {$self->{bin_count}[$i] = mean(@mov_arr);}
		elsif($fun_str eq 'median') {$self->{bin_count}[$i] = median(@mov_arr);}
	}
}

# Create sliding window count vector.
sub create_window_count{
	my $self = shift;
	if(@_ < 1) {warn "Number of bins per window not specified. Exit.\n"; return;}
	my $bins_per_window = shift;
	if(@{$self->{bin_count}} < $bins_per_window){
		warn "Bin count vector not long enough to calculate window count.\n";
		return;
	}
	# Initialize window count vector.
	$self->{window_num} = $self->{bin_num} - $bins_per_window + 1;
	$#{$self->{window_count}} = $self->{window_num} - 1;
	# $mov_window_count is the moving window count.
	my $mov_window_count = 0;
	for my $i(0..$bins_per_window-1) {$mov_window_count += $self->{bin_count}[$i];}
	$self->{window_count}[0] = $mov_window_count;
	# $i is the index to current window.
	for(my $i=1; $i < $self->{window_num}; $i++){
		$mov_window_count = $mov_window_count - $self->{bin_count}[$i-1] + $self->{bin_count}[$i-1+$bins_per_window];
		$self->{window_count}[$i] = $mov_window_count;
	}
}

# Create sliding window count vector with direction.
sub create_window_count_direction{
	my $self = shift;
	if(@_ < 1) {warn "Number of bins per window not specified. Exit.\n"; return;}
	my $bins_per_window = shift;
	if(@{$self->{bin_count_F}} < $bins_per_window or @{$self->{bin_count_R}} < $bins_per_window){
		warn "Bin count vector not long enough to calculate window count.\n";
		return;
	}
	# Initialize window count vector.
	$self->{window_num} = $self->{bin_num} - $bins_per_window + 1;
	$#{$self->{window_count_F}} = $self->{window_num} - 1;
	$#{$self->{window_count_R}} = $self->{window_num} - 1;
	# $mov_window_count is the moving window count.
	my $mov_window_count_F = 0;
	my $mov_window_count_R = 0;
	for my $i(0..$bins_per_window-1) {$mov_window_count_F += $self->{bin_count_F}[$i];}
	for my $i(0..$bins_per_window-1) {$mov_window_count_R += $self->{bin_count_R}[$i];}
	$self->{window_count_F}[0] = $mov_window_count_F;
	$self->{window_count_R}[0] = $mov_window_count_R;
	# $i is the index to current window.
	for(my $i=1; $i < $self->{window_num}; $i++){
		$mov_window_count_F = $mov_window_count_F - $self->{bin_count_F}[$i-1] + $self->{bin_count_F}[$i-1+$bins_per_window];
		$mov_window_count_R = $mov_window_count_R - $self->{bin_count_R}[$i-1] + $self->{bin_count_R}[$i-1+$bins_per_window];
		$self->{window_count_F}[$i] = $mov_window_count_F;
		$self->{window_count_R}[$i] = $mov_window_count_R;
	}
}

# Delete chromosome file but keep memory data.
sub del_file{
	my $self = shift;
	if(defined $self->{file_name}){
		unlink $self->{file_name};
	}
}

1;

__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

MyShortRead::ChromBed - A perl class to deal with single chromosome BED file.

=head1 SYNOPSIS

  use MyShortRead::ChromBed;
  my $chrBed = new MyShortRead::ChromBed;

=head1 DESCRIPTION

  This class is used to represent and perform operations on single chromosome BED
  file. After the object has been created, it can bin the genome and return a count
  vector. Then it can also calculate sliding window counts based on the bin count
  vector. Before the object dies, you should call member function to delete the BED file.

=head2 EXPORT

  This is an object-oriented module.


=head1 SEE ALSO

  MyShortRead::MyShortRead
  MyShortRead::SRBed

Mailing list: https://groups.google.com/forum/#!forum/diffreps-discuss

Web site: https://code.google.com/p/diffreps/

=head1 AUTHOR

Li Shen, E<lt>shenli.sam@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010-2013 by Li Shen

diffReps goes under GNU GPL v3: http://www.gnu.org/licenses/gpl.html


=cut
