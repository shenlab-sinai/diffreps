# Package to define genome sizes, statistics, etc.
package PJ::Genome;

use strict;
use 5.006;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use MyBioinfo::Common ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
	get_gz get_effgz get_chrsz
) ] );

our @EXPORT_OK = @{ $EXPORT_TAGS{'all'} };

our @EXPORT = ();

our $VERSION = '0.10';

# Genome sizes, chromosome sizes...
our %genome = (
	'mm9' => {
		'gsize' => 2.65e9,
		'effsize' => 1.87e9,
		'chrsize' => {
			chr1 => 197195431,
			chr2 => 181748086,
			chr3 => 159599782,
			chr4 => 155630119,
			chr5 => 152537258,
			chr6 => 149517036,
			chr7 => 152524552,
			chr8 => 131738870,
			chr9 => 124076171,
			chr10 => 129993254,
			chr11 => 121843855,
			chr12 => 121257529,
			chr13 => 120284311,
			chr14 => 125194863,
			chr15 => 103494973,
			chr16 => 98319149,
			chr17 => 95272650,
			chr18 => 90772030,
			chr19 => 61342429,
			chrX => 166650295,
			chrY => 15902554,
			chrM => 16299
		}
	},
	'hg19' => {
		'gsize' => 3.2e9,
		'effsize' => 2.7e9,
		'chrsize' => {
			chr1 => 249250621,
			chr2 => 243199373,
			chr3 => 198022430,
			chr4 => 191154276,
			chr5 => 180915260,
			chr6 => 171115067,
			chr7 => 159138663,
			chr8 => 146364022,
			chr9 => 141213431,
			chr10 => 135534747,
			chr11 => 135006516,
			chr12 => 133851895,
			chr13 => 115169878,
			chr14 => 107349540,
			chr15 => 102531392,
			chr16 => 90354753,
			chr17 => 81195210,
			chr18 => 78077248,
			chr19 => 59128983,
			chr20 => 63025520,
			chr21 => 48129895,
			chr22 => 51304566,
			chrX => 155270560,
			chrY => 59373566,
			chrM => 16571
		}
	},
	'rn4' => {
		'gsize' => 2.8e9,
		'effsize' => 1.96e9,
		'chrsize' => {
			chr1 => 267910886,
			chr2 => 258207540,
			chr3 => 171063335,
			chr4 => 187126005,
			chr5 => 173096209,
			chr6 => 147636619,
			chr7 => 143002779,
			chr8 => 129041809,
			chr9 => 113440463,
			chr10 => 110718848,
			chr11 => 87759784,
			chr12 => 46782294,
			chr13 => 111154910,
			chr14 => 112194335,
			chr15 => 109758846,
			chr16 => 90238779,
			chr17 => 97296363,
			chr18 => 87265094,
			chr19 => 59218465,
			chr20 => 55268282,
			chrX => 160699376,
			chrM => 16300
		}
	}
);

# Get genome size.
sub get_gz{
	my $gname = shift;
	if(exists $genome{$gname}){
		return $genome{$gname}->{'gsize'};
	}else{
		warn "Unavailable genome name: $gname. Return -1.\n";
		return -1;
	}
}

# Get effective genome size.
sub get_effgz{
	my $gname = shift;
	if(exists $genome{$gname}){
		return $genome{$gname}->{'effsize'};
	}else{
		warn "Unavailable genome name: $gname. Return -1.\n";
		return -1;
	}
}

# Get chromosome sizes.
sub get_chrsz{
	my $gname = shift;
	if(exists $genome{$gname}){
		return $genome{$gname}->{'chrsize'};
	}else{
		unless($gname eq ''){
			warn "Unavailable genome name: $gname. Return empty table.\n";
		}
		return {};
	}
}


1;

__END__

=head1 NAME

PJ::Genome - Module for convenient access to genome statistics, such as chromosome lengths.

=head1 SYNOPSIS

  use PJ::Genome qw( get_gz get_effgz get_chrsz );

  my $genome_size = &get_gz($genome_name);  # genome size.
  my $eff_genome_size = &get_effgz($genome_name);  # effective genome size.
  my $ref_hash_chrom_size = &get_chrsz($genome_name);  # chromosome lengths in a reference to hash table.


=head1 DESCRIPTION

  Many bioinformatics functions need information about a genome, such as the genome size,
  chromosome lengths. This module contains these information for a few model species and
  provide PERL interfaces for easy access.

=head2 EXPORT

  None.

=head1 SEE ALSO

MyBioinfo::Common

Mailing list: https://groups.google.com/forum/#!forum/diffreps-discuss

Web site: https://code.google.com/p/diffreps/


=head1 AUTHOR

Li Shen, E<lt>shenli.sam@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010-2013 by Li Shen

diffReps goes under GNU GPL v3: http://www.gnu.org/licenses/gpl.html


=cut
