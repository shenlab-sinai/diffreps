#
# My Perl library to deal with short read sequencing data.
# created by Li Shen, Nov 2009.
#
# strGCPercent added March 2010 - Li Shen
# bin_genome_count,read_chrom_len,sep_chrom_bed,del_chrfile added April 2010 - Li Shen
# bin_genome_count_direction added June 2010 - Li Shen
#
#

use 5.006;
use strict;
use warnings;
package MyShortRead::MyShortRead;


use POSIX;
use Time::HiRes qw(time);
use File::Temp qw(tempfile);

require Exporter;
our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use MyShortRead::MyShortRead ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(

) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(compare2 overlap2 strGCPercent bin_genome_count bin_genome_count_direction read_chrlen_tbl read_chrlen_ordered sep_chrom_bed del_chrfile order_chr);

our $VERSION = '1.00';


# Preloaded methods go here.

# Comparing two genomic intervals according to start or end position.
sub compare2{
  my($a,$b) = @_;
  my $cmpby = 'start';
  if(@_ > 2) {$cmpby = $_[2];}

  $a->{chrom} =~ /^(chr)?(\w+)$/;
  my $chr1 = $2;
  $b->{chrom} =~ /^(chr)?(\w+)$/;
  my $chr2 = $2;

  # $chr1 = 100 if uc $chr1 eq 'X';
  # $chr1 = 101 if uc $chr1 eq 'Y';
  # $chr1 = 102 if uc $chr1 eq 'M';
  # $chr1 = 199 if uc $chr1 eq 'Z';
  # $chr2 = 100 if uc $chr2 eq 'X';
  # $chr2 = 101 if uc $chr2 eq 'Y';
  # $chr2 = 102 if uc $chr2 eq 'M';
  # $chr2 = 199 if uc $chr2 eq 'Z';

  # return $chr1 <=> $chr2 unless $chr1 == $chr2;

  return $chr1 cmp $chr2 unless $chr1 eq $chr2;

  if($cmpby eq 'end') {return $a->{'end'} <=> $b->{'end'};}
  else {return $a->{'start'} <=> $b->{'start'};}
}

# Order chromosome names.
sub order_chr{
  my($a,$b) = @_;
  my($chr1, $chr2);
  # Find chromosome names.
  if(exists $a->{'chr'}) {$chr1 = $a->{'chr'};}
  elsif(exists $a->{chrom}) {$chr1 = $a->{chrom};}
  else {warn "Cannot find chromosome name! Order is not determined.\n"; return 0;}
  if(exists $b->{'chr'}) {$chr2 = $b->{'chr'};}
  elsif(exists $b->{chrom}) {$chr2 = $b->{chrom};}
  else {warn "Cannot find chromosome name! Order is not determined.\n"; return 0;}

  return $chr1 cmp $chr2;

  # # Get rid of 'chr' if exists.
  # $chr1 =~ /^(chr)?(\w+)$/; $chr1 = $2;
  # $chr2 =~ /^(chr)?(\w+)$/; $chr2 = $2;
  # # Convert letters to large numbers.
  # $chr1 = 100 if uc $chr1 eq 'X';
  # $chr1 = 101 if uc $chr1 eq 'Y';
  # $chr1 = 102 if uc $chr1 eq 'M';
  # $chr1 = 199 if uc $chr1 eq 'Z';
  # $chr2 = 100 if uc $chr2 eq 'X';
  # $chr2 = 101 if uc $chr2 eq 'Y';
  # $chr2 = 102 if uc $chr2 eq 'M';
  # $chr2 = 199 if uc $chr2 eq 'Z';

  # return $chr1 <=> $chr2;
}

# Find the interval list with the lowest start or end position.
# return the index of the lowest list and set the position through reference.
sub find_lowest_list{
  my($rlists, $rpos) = @_;
  my $cmpby = 'start';
  if(@_ > 2) {$cmpby = $_[2];}
  my $ll = 0; # initialize lowest list.
  my $lowest_intrv = $rlists->[$ll]{cur_intrv}; # intrv with the lowest end position.
  for my $i(1..$#{$rlists}){
    my $cur_intrv = $rlists->[$i]{cur_intrv};
    if(compare2($cur_intrv, $lowest_intrv, $cmpby) < 0){
      $lowest_intrv = $cur_intrv;
      $ll = $i;
    }
  }
  # return the lowest position through reference.
  ${$rpos} = $cmpby eq 'end'? $lowest_intrv->{'end'} : $lowest_intrv->{'start'};
  return $ll;
}

# Judge whether two intervals overlap.
sub overlap2{
  my($refintrv1, $refintrv2) = @_;
  if($refintrv1->{chrom} eq $refintrv2->{chrom}){
    if($refintrv1->{'start'} <= $refintrv2->{'end'} and
       $refintrv1->{'end'} >= $refintrv2->{'start'}){
      return 1;
    }
  }
  return 0;
}

# Calculate the GC percentage given a DNA sequence in string representation.
sub strGCPercent{
  my $dnaSeq = shift;
  my $c = ($dnaSeq =~ s/[gc]/X/gi); # Count the # of G or C without regarding to case.
  my $n = ($dnaSeq =~ s/n/Y/gi);	# count the # of 'N's.
  return ($c+$n/2) / length($dnaSeq);
}

# Given a BED file, bin the genome and count the reads in bins.
sub bin_genome_count{
	# BED file,bin size,chromosome length,fragment size,reference to bin vector.
	my($bedfile,$binsize,$chrlen,$fragsize,$r_binVec) = @_;

    $#{$r_binVec} = ceil($chrlen/$binsize)-1;   # Calculate bin vector size.
    foreach my $i(0..$#{$r_binVec}) {$r_binVec->[$i] = 0;}  # Initialize bin vector to zeros.

	unless(open HBED, "<", $bedfile) {
        warn "Open bed file $bedfile error: $!\n";
        return;
    }
	my $line = 0;
	while(<HBED>){
		$line++;
		chomp;
		my @cells = split /\t/;
		if(@cells < 6){
			warn "Format error at $bedfile, line:$line. Skip!\n";
			next;
		}
		my($chrom,$start,$end,$name,$score,$strand) = @cells;
		# Shift read location by fragment size.
		my $loc;
		if($strand eq '+') {$loc = $start + $fragsize/2;}
		else {$loc = $end - $fragsize/2;}
		# Increment read count in bin.
		# $loc--;	# Decrement location to correctly assign the read to bin vector.(NO! BED is 0-based.)
		$r_binVec->[floor($loc/$binsize)]++;
	}
	close HBED;
}

# Bin the genome and count the reads with direction.
sub bin_genome_count_direction{
	# BED file,bin size,chromosome length,fragment size,reference to bin vector.
	my($bedfile,$binsize,$chrlen,$fragsize,$r_binVec_F,$r_binVec_R) = @_;

    $#{$r_binVec_F} = ceil($chrlen/$binsize)-1; # Calculate bin vector size.\
    $#{$r_binVec_R} = ceil($chrlen/$binsize)-1; # Calculate bin vector size.\
    foreach my $i(0..$#{$r_binVec_F}) {$r_binVec_F->[$i] = 0;}  # Initialize bin vector to zeros.
    foreach my $i(0..$#{$r_binVec_R}) {$r_binVec_R->[$i] = 0;}  # Initialize bin vector to zeros.

	unless(open HBED, "<", $bedfile) {
        warn "Open bed file $bedfile error: $!\n";
        return;
    }
	my $line = 0;
	while(<HBED>){
		$line++;
		chomp;
		my @cells = split;
		if(@cells < 6){
			warn "Format error at $bedfile, line:$line. Skip!\n";
			next;
		}
		my($chrom,$start,$end,$name,$score,$strand) = @cells;
		# Shift read location by fragment size and then increment the read count at the correct location.
		my $loc;
		if($strand eq '+') {
			$loc = $start + $fragsize/2;
			$r_binVec_F->[floor($loc/$binsize)]++;
		}
		else {
			$loc = $end - $fragsize/2;
			$r_binVec_R->[floor($loc/$binsize)]++;
		}
	}
	close HBED;
}

# Read chromosome name and length into a hash table.
sub read_chrlen_tbl{
	my $chrfile = shift;
	open HCHR, "<", $chrfile or die "Cannot open chromosome length file:$!\n";
	my $line = 0;
	my %chrlen;	# hash table for chromosome length.
	while(<HCHR>){
		$line++;
		chomp;
		next if $_ eq '';	# skip empty line.
		my @cells = split;
		warn "Insufficient chromosome information at $chrfile, line:$line. Skip!\n" if @cells < 2;
		my($chrom,$len) = @cells;
		if($chrom =~ /^\S+$/ and $len =~ /^[0-9]+$/){
			$chrlen{$chrom} = $len;
		}
		else {warn "Format error at $chrfile, line:$line. Skip!\n"}
	}
	return %chrlen;
}

# Read chromosome name and length into an array and order them by name.
sub read_chrlen_ordered{
	my $chrfile = shift;
	open HCHR, "<", $chrfile or die "Cannot open chromosome length file:$!\n";
	my $line = 0;
	my @chrlen;	# array for chromosome length.
	while(<HCHR>){
		$line++;
		chomp;
		next if $_ eq '';	# skip empty line.
		my @cells = split;
		warn "Insufficient chromosome information at $chrfile, line:$line. Skip!\n" if @cells < 2;
		my($chrom,$len) = @cells;
		if($chrom =~ /^\S+$/ and $len =~ /^[0-9]+$/){
			my %rec = (
				chrom	=> $chrom,
				len		=> $len);
			push @chrlen, \%rec;
		}
		else {warn "Format error at $chrfile, line:$line. Skip!\n"}
	}
	return sort {order_chr($a,$b)} @chrlen;	# return ordered list.
}


# Separate a BED file into small ones by chromosome names.
sub sep_chrom_bed{
	my($bedfile,@chrname) = @_;
	my $timestamp = time;
	my $timehead = ".$timestamp.";	# timestamp file name head.
	my %hsh_chrfile;	# hash of chromosome file handles.
	my %hsh_chrname;	# hash of chromosome file names.
	# Create chromosome files for writing.
	foreach my $chrom(@chrname){
		my ($fh, $chrfile) = tempfile( "diffreps$chrom.bed.XXXXX", TMPDIR => 1 );
		$hsh_chrfile{$chrom} = $fh;
		$hsh_chrname{$chrom} = $chrfile;
	}
	open HBED, "<", $bedfile or die "Cannot open $bedfile:$!\n";
	my $line = 0;	# line number in BED file.
	my $nread = 0;	# actual number of reads splitted.
	# Go through bed file and assign each line to corresponding file.
	while(<HBED>){
		$line++;
		chomp;
		next if $_ eq '';
		my @cells = split;
		if(@cells < 6){
			warn "Insufficient tag information at $bedfile, line:$line. Skip!\n";
			next;
		}
		my($chrom,$start,$end,$name,$score,$strand) = @cells;
		if(exists $hsh_chrfile{$chrom}){
			$nread++;
			print {$hsh_chrfile{$chrom}} $_ . "\n";
		}
		else{
			warn "Unspecified chromosome name at $bedfile, line:$line.\nCheck your chromosome length file. Skip!\n";
		}
	}
	# Close all files.
	close HBED;
	foreach my $hdl(values %hsh_chrfile) {close $hdl;}

	# Return chromosome file name vector.
	return ($nread, \%hsh_chrname);
}

# Delete all separated chromosome files.
sub del_chrfile{
	foreach my $file(@_) {
		unlink $file or warn "Cannot delete $file: $!\n";
	}
}

1;

__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

MyShortRead::MyShortRead - My Perl library to deal with nextgen short read data.

=head1 SYNOPSIS

  use MyShortRead::MyShortRead;

  my $compare_res = &compare2($interval_a, $interval_b, $compare_by);
  my $boolean = &overlap2(\@interval1, \@interval2);
  my $GC_percent = &strGCPercent($DNA_seq);
  &bin_genome_count($bedfile, $bin_size, $chrom_len, $fragsize, \@bin_vector);
  &bin_genome_count_direction($bedfile, $bin_size, $chrom_len, $fragsize, \@bin_vector);
  my %chrom_tbl = &read_chrlen_tbl($chrom_filename);
  my @chrom_list = &read_chrlen_ordered($chrom_filename);
  my ($n_reads, $ref_hash_chrom_sepFileName) = &sep_chrom_bed($bedfile, @chrom_name);
  &del_chrfile(@chrom_sepFiles);
  my $compare_res = &order_chr($chrom_name_a, $chrom_name_b);

=head1 DESCRIPTION

  This module contains a few functions that I created to perform some operations on
  short read data from next-generation sequencing machines. They include separate a
  BED file according to chromosomes; bin genome and count the number of reads; compare
  two genomic intervals and determine which one comes first(on the left side).

=head2 EXPORT

  compare2
  overlap2
  strGCPercent
  bin_genome_count
  bin_genome_count_direction
  read_chrlen_tbl
  read_chrlen_ordered
  sep_chrom_bed
  del_chrfile
  order_chr



=head1 SEE ALSO

MyShortRead::ChromBed
MyShortRead::SRBed


=head1 AUTHOR

Li Shen, E<lt>shenli.sam@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010-2013 by Li Shen

diffReps goes under GNU GPL v3: http://www.gnu.org/licenses/gpl.html


=cut

