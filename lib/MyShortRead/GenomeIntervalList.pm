#!/usr/bin/perl -w
# GenomeIntervalList - a class for manipulation of genomic interval list.
# by Li Shen April 2010.
#

use strict;
use 5.006;
package MyShortRead::GenomeIntervalList;
our $VERSION = '0.01';


# Class constructor.
sub new{
  my $class = shift;
  my $self = {};
  $self->{file} = undef;
  $self->{list} = [];
  $self->{cur_intrv} = undef;
  $self->{cur_pos} = undef;
  $self->{sum_info} = undef;
  bless($self, $class);
  
  # user supplied file name and structure.  
  # do nothing if user does not supply enough information.
  if(@_ >= 6) {$self->readfile(@_);}
 
  return $self; 
}

# Sum info in the list.
sub sum_info{
  my $self = shift;
  $self->{sum_info} = 0;
  foreach my $ri(@{$self->{list}}){
    $self->{sum_info} += $ri->{'info'} if $ri->{'chr'} ne 'Z'; # avoid infinite interval.
  }
}

# Read genome interval list from file.
sub readfile{
  my $self = shift;
  # File name, header, column index for chromosome, start, end, info.
  my($intrvfile,$header,$chr_c,$start_c,$end_c,$info_c) = @_;
  my @onelist; # create an empty array for a list.
  open HINTRV, '<', $intrvfile or die "Open interval file $intrvfile error: $!\n";
  <HINTRV> if $header; # Skip header.
  my $line = 0;
  while(<HINTRV>){
    $line++;
    chomp;
    my @cells = split;
    my %arec; # create an empty hash for a record.
    die "chromosome position format unrecognized at line: $line.\n" unless $cells[$chr_c] =~ /^(chr|c)?(\d+|[xym])$/i;
    next if uc $2 eq 'M'; # Ignore chrM.
    $arec{'chr'} = $2;
    $arec{'start'} = $cells[$start_c];
    $arec{'end'} = $cells[$end_c];
    $arec{'info'} = $cells[$info_c];
    push @onelist, \%arec;
  }
  close HINTRV;
  $self->{list} = \@onelist; # store the list as object member.
  $self->{file} = \$intrvfile; # book-keep file name.
  $self->{cur_pos} = 0; # reset current position.
  $self->{cur_intrv} = $self->{list}[0];
}

# Sort interval list if needed. Add an infinite interval to the list end.
sub sort_addinf{
  my $self = shift;
  my $needsort = 0; # default = no sort.
  if(@_ > 0) {$needsort = shift;} # Boolean indicates sorting.
  if($needsort){
    @{$self->{list}} = sort cmploc @{$self->{list}};
    $self->reset();
  }
  # Check if the last interval is already an infinite.
  my $last_intrv = $self->{list}[-1]; 
  if($last_intrv->{'chr'} ne 'Z'){
    # An infinite interval never overlaps with any interval(incl. itself).
    my $infintrv = {}; 
    $infintrv->{'chr'} = 'Z';
    $infintrv->{'start'} = 1;
    $infintrv->{'end'} = -1; # let end < start.
    push @{$self->{list}}, $infintrv;
  }
}

# Compare function for sorting genomic intervals according to start position. 
# Use 100 for chrX and 101 for chrY for ease of comparison.
# Use 199 for chrZ so the infinite interval is always ranked last.
sub cmploc{
  my $chr1 = $a->{'chr'};
  my $chr2 = $b->{'chr'};
  $chr1 = 100 if uc $chr1 eq 'X';
  $chr1 = 101 if uc $chr1 eq 'Y';
  $chr1 = 102 if uc $chr1 eq 'M';
  $chr1 = 199 if uc $chr1 eq 'Z';
  $chr2 = 100 if uc $chr2 eq 'X';
  $chr2 = 101 if uc $chr2 eq 'Y';
  $chr2 = 102 if uc $chr2 eq 'M';
  $chr2 = 199 if uc $chr2 eq 'Z';

  return $chr1 <=> $chr2 unless $chr1 == $chr2;
  return $a->{'start'} <=> $b->{'start'};
}

# Advance the list by one position. Do nothing if already at the end.
sub advance{
  my $self = shift;
  if($self->{cur_intrv}{'chr'} ne 'Z'){
    $self->{cur_pos}++;
    $self->{cur_intrv} = $self->{list}[$self->{cur_pos}];
  }
} 

# Reset the current position to the start.
sub reset{
  my $self = shift;
  $self->{cur_pos} = 0;
  $self->{cur_intrv} = $self->{list}[0];
}

# Tell whether the current position is at the end.
sub is_end{
  my $self = shift;
  return 1 if $self->{cur_intrv}{'chr'} eq 'Z';
  return 0;
}

1;
