# SRBed - A class for manipulation of a short read dataset in BED file.
# Bins - are non-overlapping windows that span the whole genome.
# Windows - are sliding windows with fixed size. A sliding step equals the bin size.
# Sliding window: when implementing a sliding window procedure, the whole genome is
# first separated into non-overlapping bins of the size of the sliding step. Then a 
# window count is calculated for a continuous block of bins.
# window_size must be multiple times of bin_size, i.e. window_size % bin_size = 0.


use strict;
use 5.006;
package MyShortRead::SRBed;
our $VERSION = '1.00';

use MyShortRead::MyShortRead qw(sep_chrom_bed);
use MyShortRead::ChromBed;

# Class constructor.
sub new{
	my $class = shift;
	my $self = {};
	$self->{file_name} = undef;
	$self->{bin_size} = 1000;
	$self->{window_size} = 1000;
	$self->{frag_size} = 100;
	$self->{chr_list} = {};	# hash table: chrom -> ChromBed
	bless($self, $class);

	# User supplied bed file name.
	if(@_ >= 1) {$self->{file_name} = $_[0];}
	if(@_ >= 2) {$self->{bin_size} = $_[1];}
	if(@_ >= 3){$self->{frag_size} = $_[2];}

	return $self;
}

sub set_file_name{
	my $self = shift;
	if(@_ >= 1) {$self->{file_name} = $_[0];}
}

sub set_bin_size{
	my $self = shift;
	if(@_ >= 1) {$self->{bin_size} = $_[0];}
}

sub set_window_size{
	my $self = shift;
	if(@_ >= 1) {$self->{window_size} = $_[0];}
}

sub set_frag_size{
	my $self = shift;
	if(@_ >= 1) {$self->{frag_size} = $_[0];}
}

# Create hash table of chrom name to ChromBed object given a chrlen hash table. 
sub create_chr_list{
	my $self = shift;
	my $rtbl = $_[0];	# chrlen hash table.
	# Separate the whole bed file into trunks of chromosomes.
	my($nread, $hsh_chrfile) = sep_chrom_bed($self->{file_name}, keys %{$rtbl});
	# For each chromosome, create a ChromBed object.
	while(my($chrom,$chrfile) = each %{$hsh_chrfile}){
		my $chrlen = $rtbl->{$chrom};
		my $chrbed = new MyShortRead::ChromBed($chrom,$chrlen,$chrfile);
		# $chrbed->do_bin_count($self->{bin_size},$self->{frag_size});
		$self->{chr_list}{$chrom} = $chrbed;	# add the new obj to hash table.
	}
	return $nread;
}

# Do bin count for one chromosome.
sub chr_bin_count{
	my $self = shift;
	my $chrom = shift;
	if(exists $self->{chr_list}{$chrom}){
		$self->{chr_list}{$chrom}->do_bin_count($self->{bin_size},$self->{frag_size});
	}
	else {warn "Chromosome name not exists! Bin count not done.\n";}
}

# Do bin count for one chromosome with direction.
sub chr_bin_count_direction{
	my $self = shift;
	my $chrom = shift;
	if(exists $self->{chr_list}{$chrom}){
		$self->{chr_list}{$chrom}->do_bin_count_direction($self->{bin_size},$self->{frag_size});
	}
	else {warn "Chromosome name not exists! Bin count not done.\n";}
}

# Smooth bin count signals for one chromosome.
sub chr_bin_smooth{
	my $self = shift;
	if(@_ < 3){
		warn "Not enough arguments: chromosome name, smooth width and function string must be specified.\n";
		return;
	}
	my($chrom,$smooth_width,$fun_str) = @_;
	if(exists $self->{chr_list}{$chrom}){
		$self->{chr_list}{$chrom}->smooth_bin_count($smooth_width,$fun_str);
	}
	else {warn "Chromosome name not exists! Bin smooth not done.\n";}
}

# Calculate sliding window count for one chromosome.
sub chr_win_count{
	my $self = shift;
	my $chrom = shift;
	my $bins_per_window = $self->{window_size} / $self->{bin_size};
	if(exists $self->{chr_list}{$chrom}){
		if(defined $self->{chr_list}{$chrom}->{bin_num}){
			$self->{chr_list}{$chrom}->create_window_count($bins_per_window);
		}
		else {warn "Bin count vector must be created before window count.\n";}
	}
	else {warn "Chromosome name not exists! Window count not done.\n";}
}

# Calculate sliding window count for one chromosome with direction.
sub chr_win_count_direction{
	my $self = shift;
	my $chrom = shift;
	my $bins_per_window = $self->{window_size} / $self->{bin_size};
	if(exists $self->{chr_list}{$chrom}){
		if(defined $self->{chr_list}{$chrom}->{bin_num}){
			$self->{chr_list}{$chrom}->create_window_count_direction($bins_per_window);
		}
		else {warn "Bin count vector must be created before window count.\n";}
	}
	else {warn "Chromosome name not exists! Window count not done.\n";}
}

# Delete bin count vector for specified chromosome.
sub del_bin_count{
	my $self = shift;
	my $chrom = shift;
	if(exists $self->{chr_list}{$chrom}){
		$self->{chr_list}{$chrom}->{bin_count} = [];
		$self->{chr_list}{$chrom}->{bin_count_F} = [];
		$self->{chr_list}{$chrom}->{bin_count_R} = [];
	}
}

# Delete window count vector for specified chromosome.
sub del_win_count{
	my $self = shift;
	my $chrom = shift;
	if(exists $self->{chr_list}{$chrom}){
		$self->{chr_list}{$chrom}->{window_count} = [];
		$self->{chr_list}{$chrom}->{window_count_F} = [];
		$self->{chr_list}{$chrom}->{window_count_R} = [];
	}
}

# Return read count given chrom name and bin position.
sub get_bin_count{
	my $self = shift;
	my($chrom,$bin_pos) = @_;
	if(exists $self->{chr_list}{$chrom}){
		my $r_chrbed = $self->{chr_list}{$chrom};
		if($bin_pos < $r_chrbed->{bin_num}){
			return $r_chrbed->{bin_count}->[$bin_pos];
		}
		else{
			my $bin_num = $r_chrbed->{bin_num};
			warn "Chrom: $chrom, Bin position:$bin_pos exceeds boundary:$bin_num! Bin count zero returned.\n";
			return 0;
		}
	}
	else{
		warn "Cannot find chrom name! Bin count zero returned.\n";
		return 0;
	}
}

# Return read count with direction given chrom name and bin position.
sub get_bin_count_direction{
	my $self = shift;
	my($chrom,$bin_pos) = @_;
	if(exists $self->{chr_list}{$chrom}){
		my $r_chrbed = $self->{chr_list}{$chrom};
		if($bin_pos < $r_chrbed->{bin_num}){
			return ($r_chrbed->{bin_count_F}->[$bin_pos], $r_chrbed->{bin_count_R}->[$bin_pos]);
		}
		else{
			my $bin_num = $r_chrbed->{bin_num};
			warn "Chrom: $chrom, Bin position:$bin_pos exceeds boundary:$bin_num! Bin count zero returned.\n";
			return 0;
		}
	}
	else{
		warn "Cannot find chrom name! Bin count zero returned.\n";
		return 0;
	}
}

# Return read count given chrom name and window position.
sub get_win_count{
	my $self = shift;
	my($chrom,$win_pos) = @_;
	if(exists $self->{chr_list}{$chrom}){
		my $r_chrbed = $self->{chr_list}{$chrom};
		if($win_pos < $r_chrbed->{window_num}){
			return $r_chrbed->{window_count}[$win_pos];
		}
		else{
			warn "Window position exceeds boundary! Window count zero returned.\n";
			return 0;
		}
	}
	else{
		warn "Cannot find chrom name! Window count zero returned.\n";
		return 0;
	}
}

# Return read count given chrom name and window position.
sub get_win_count_direction{
	my $self = shift;
	my($chrom,$win_pos) = @_;
	if(exists $self->{chr_list}{$chrom}){
		my $r_chrbed = $self->{chr_list}{$chrom};
		if($win_pos < $r_chrbed->{window_num}){
			return ($r_chrbed->{window_count_F}[$win_pos], $r_chrbed->{window_count_R}[$win_pos]);
		}
		else{
			warn "Window position exceeds boundary! Window count zero returned.\n";
			return 0;
		}
	}
	else{
		warn "Cannot find chrom name! Window count zero returned.\n";
		return 0;
	}
}

# Return read count given a genomic region's coordinates.
sub get_region_count{
	my $self = shift;
	my($chrom,$start,$end) = @_;
	if(($start-1) % $self->{bin_size} != 0 or $end % $self->{bin_size} != 0){
		warn "start or end coordinate not on bin boundary. round to the nearest bin!\n";
	}
	my $bin_start = int(($start-1)/$self->{bin_size});
	my $bin_end = int($end/$self->{bin_size}) - 1;
	my $count = 0;
	if($bin_start > $bin_end) {warn "Genomic region contains zero bin! Count zero returned.\n";}
	else {for my $i($bin_start..$bin_end) {$count += $self->get_bin_count($chrom,$i);} }
	return $count;
}

# Return read count with direction given a genomic region's coordinates.
sub get_region_count_direction{
	my $self = shift;
	my($chrom,$start,$end) = @_;
	if(($start-1) % $self->{bin_size} != 0 or $end % $self->{bin_size} != 0){
		warn "start or end coordinate not on bin boundary. round to the nearest bin!\n";
	}
	my $bin_start = int(($start-1)/$self->{bin_size});
	my $bin_end = int($end/$self->{bin_size}) - 1;
	my @count = (0,0);	# (Forward, Reverse).
	if($bin_start > $bin_end) {warn "Genomic region contains zero bin! Count zero returned.\n";}
	else {
		for my $i($bin_start..$bin_end) {
			my @bin_count = $self->get_bin_count_direction($chrom,$i);
			$count[0] += $bin_count[0];
			$count[1] += $bin_count[1];
		} 
	}
	return @count;
}

# Delete all associated chromosome files but keep memory data.
sub del_all_chrfiles{
	my $self = shift;
	return if keys %{$self->{chr_list}} == 0;
	while(my($chrom,$chrbed) = each %{$self->{chr_list}}){
		$chrbed->del_file();
	}
}

1;


__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

MyShortRead::SRBed - A perl class to deal with a short read dataset in BED file.

=head1 SYNOPSIS

  use MyShortRead::SRBed;
  my $srBed = new MyShortRead::SRBed;

=head1 DESCRIPTION

  This class is used to extract information from a short read BED file. It first
  separates the whole file into small pieces of chromosome files and then do bin
  and sliding window count. Before the object exists, you should call its member
  function to delete all temporary chromosome files.


=head2 EXPORT

  This is an object-oriented module.


=head1 SEE ALSO

  MyShortRead::MyShortRead
  MyShortRead::ChromBed


=head1 AUTHOR

Li Shen, E<lt>shenli.sam@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010-2013 by Li Shen

diffReps goes under GNU GPL v3: http://www.gnu.org/licenses/gpl.html


=cut



