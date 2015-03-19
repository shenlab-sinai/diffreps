#!/usr/bin/perl -w
# A class to read and parse Illumina Fastq files.

use strict;
use 5.006;
package MyShortRead::Fastq;
our $VERSION = '0.10';

# Constructor.
sub new{
	my($class) = @_;
	my $self = {
		file => undef,	# fastq file name.
		h => undef,	# file handle.
		ln => undef	# current line number in fastq file.
	};
	bless $self, $class;
	if(@_ > 1){
		$self->openFile($_[1]);
	}
	return $self;
}

# Is end of Fastq file?
sub isEnd{
	my($self) = @_;
	return eof $self->{h};
}

# Open a Fastq file.
sub openFile{
	my($self,$file) = @_;
	if($file eq '-') {$self->{h} = *STDIN;}
	else {open $self->{h}, '<', $file or die "Open Fastq file error: $!\n";}
	$self->{file} = $file;
	$self->{ln} = 0;
	return $self->{h};
}

# Read one short read sequence and quality. Return data in hash table.
sub readSeq{
	my($self) = @_;
	my $c = 0;
	my @sr;
	while(readline $self->{h}){
		chomp;
		push @sr, $_;
		$self->{ln}++;
		last if $c++ >= 3;
	}
	if($c == 4){
		if($sr[0] =~ /^@(.+):([0-9]):([0-9]+):([0-9]+):([0-9]+)#(.+)\/([1-2])$/){
			my($instr,$lane,$tile,$xcoo,$ycoo,$index,$pair) = ($1,$2,$3,$4,$5,$6,$7);
			my $sq = {
				instrument => $instr,	# instrument name.
				lane => $lane,	# lane number.
				tile => $tile,	# tile numer.
				xcoor => $xcoo,	# X-coordinate.
				ycoor => $ycoo,	# Y-coordinate.
				index => $index,	# index number.
				pair => $pair,	# read pair: 0 or 1.
				seq => $sr[1],	# read sequence.
				qual => $sr[3]	# read quality string.
			};
			return $sq;
		}
		else {return {};}
	}
	else{return {};}
}

# Calculate average Phred score for a SAM/Sanger quality string.
sub avgphred_ss{
	my($self,$str) = @_;
	my $q = 0;
	while($str =~ /(.)/g) {$q += ord($1) - 33;}
	return $q / length($str);
}

# Calculate average Phred score for a Illumina quality string.
sub avgphred_is{
	my($self,$str) = @_;
	my $q = 0;
	while($str =~ /(.)/g) {$q += ord($1) - 64;}
	return $q / length($str);
}

1;
