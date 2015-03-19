package PJ::Database::Root;

use 5.006;
use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use PJ::Database::Root ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
	
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	
);

our $VERSION = '1.00';


# Preloaded methods go here.

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

PJ::Database::Root - Perl module for genomic annotations

=head1 SYNOPSIS

  This is a dummy module which contains NOTHING. It serves as a beacon for
  other Perl programs to find the path to genomic annotations.

  require PJ::Database::Root;
  my $path = $INC{'PJ/Database/Root.pm'};
  $path = s/Root.pm$//;
  $path2genome = $path . $genome_name;	# mm9, rn4, etc.


=head1 DESCRIPTION

Program such as region_analysis.pl requires genomic annotations to properly
assign functions to regions generated from peak or differential site detection
program such as diffReps. 

The genomic annotations, such as coordinates of genebody, TSS, pericentromere
are stored in text files. It is important for region_analysis.pl to find those
files and perform operations on them.

The module PJ::Database::Root is simply an empty .pm file that does nothing. 
However, it is convenient to have this module so that its path can be extracted
from %INC.


=head2 EXPORT

None by default.

=head1 SEE ALSO

PJ::Genome

Mailing list: https://groups.google.com/forum/#!forum/diffreps-discuss

Web site: https://code.google.com/p/diffreps/


=head1 AUTHOR

Li Shen, E<lt>shenli.sam@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010-2013 by Li Shen

diffReps goes under GNU GPL v3: http://www.gnu.org/licenses/gpl.html



=cut
