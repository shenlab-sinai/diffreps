#!/usr/bin/perl -w
#
# sam2bed.pl - Convert SAM format to BED6 format.
# assume all chromosome names are in the format of
# '1','2',...,'X','Y','M'. All the others are considered as non-chromosome.
# 
# History:
#
# Added an option to output only uniquely mapped reads. 07/11/2010.
# 

# Program arguments.
use strict;
if(@ARGV < 2){
	print "Usage: $0 input.sam output.bed [-u]\n\n";
	print "-u    Option to keep only uniquely aligned reads.\n";
	print "      For reads with multiple aligned locations, only the top hit will be reported.\n";
	exit 0;
}
if($ARGV[0] eq '-') {*HSAM=*STDIN;}
else {open HSAM, "<", $ARGV[0] or die "Open input file error: $!\n";}
if($ARGV[1] eq '-') {*HBED=*STDOUT;}
else {open HBED, ">", $ARGV[1] or die "Open output file error: $!\n";}

my $utag = 0;	# unique mapping tag.
$utag = 1 if defined $ARGV[2] and $ARGV[2] eq '-u';	


# Process the rest of the SAM.
while(<HSAM>){
	next if /^\@/;	# skip SAM header lines.
	next if /^\s*\#/;	# skip comments.
	chomp;
	my($chrom,$pos,$len1,$len2,$nskip,$strand,$uniq) = parseln($_);
	if($chrom eq '') { die "Encountered lines with incomplete information. Exit.\n";}
	elsif($chrom eq 'BBB') { next;}	# skip reads mapped to non-standard chromosomes.
	if(($utag and $uniq) or !$utag){
		# pos is already 0-based.
		printf HBED "%s\t%d\t%d\t%s\t%d\t%s\n", $chrom,$pos,$pos+$len1,'N','0',$strand;
		if($nskip > 0){	# 2nd half of the spliced read.
			printf HBED "%s\t%d\t%d\t%s\t%d\t%s\n", $chrom,$pos+$len1+$nskip,$pos+$len1+$nskip+$len2,'N','0',$strand;
		} 
	}
}
close HSAM;
close HBED;



###################################
# Function to parse one SAM line.
sub parseln{
	my @cells = split /\t/, $_[0];
	return('',-1,-1,'','') if @cells < 11;	# lines with incomplete information.
	#### Strand ####
	my $strand = $cells[1] & 16? '-' : '+';	# Strand. 16 is the tag for reverse strand.
	#### Chromosome name ####
	# my $chrom;
	# if($cells[2] =~ /(chr)?([0-9]+|[XYM]|MT)(\.fa)?/i){
	# 	$chrom = 'chr' . $2;	# Extract chromosome name.
	# }else{
	# 	return('BBB',-1,-1,'','') unless defined $chrom;	# reads mapped to non-chromosome locations.
	# }
	$cells[2] =~ s/^chr//;
	$cells[2] =~ s/\.fa$//;
	my $chrom = 'chr' . $cells[2];	# chromosome name.
	#### Mapped location ####
	my $pos = $cells[3]-1;	# SAM: 1-based, left-most; BED: 0-based.
	# Determine matched and skipped lengths.
	my($len1,$len2,$nskip);
	my $cigar = $cells[5];
	if($cigar =~ /^([0-9]+)M(([0-9]+)N([0-9]+)M)?$/){	# Bowtie format.
		$len1 = $1;
		if(defined $2){	# the read is spliced.
			$len2 = $4;
			$nskip = $3;
		} else{
			$len2 = 0;
			$nskip = 0;
		}
	} else{	# Ignore other formats and set to seq length.
		$len1 = length($cells[9]);
		$len2 = 0;
		$nskip = 0;
	}
	# Determine uniqueness of mapping.
	my $uniq = 1;
	if(@cells > 11){	# processing tags.
		my @tags = splice @cells, 11;
		my %tag_tbl;
		foreach my $t(@tags){
			my @e = split /\:/, $t;
			$tag_tbl{$e[0]} = {'type' => $e[1], 'value' => $e[2]};
		}
		if(exists $tag_tbl{'XT'}){	# BWA read type tag.
			if($tag_tbl{'XT'}->{'value'} ne 'U') {$uniq = 0;}
		} elsif(exists $tag_tbl{'NH'}){	# SAM tag for number of alignments: Bowtie and others.
			if($tag_tbl{'NH'}->{'value'} > 1) {$uniq = 0;}
		}
	}
	return ($chrom,$pos,$len1,$len2,$nskip,$strand,$uniq);
}

