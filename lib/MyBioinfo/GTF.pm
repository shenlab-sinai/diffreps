package MyBioinfo::GTF;

use strict;

require Exporter;
our @ISA = qw(Exporter);

our @EXPORT_OK = qw(read_gtf_tbl anno_gene_table);

our @EXPORT = qw(read_gtf_tbl anno_gene_table);

our $VERSION = '0.50';

# Read gene structures into a hash table.
sub read_gtf_tbl{
	my $gtf = shift;	# GTF file name.
	open GTF, "<$gtf" or die "Open input gtf file error:$!\n";
	my %gene_table;
	my %cur_transcript = (	# working transcript data structure.
		'chrom' => '',
		'strand' => '',
		'source' => '',
		'gid' => '',
		'gname' => '',
		'tid' => '',
		'exons' => []
	);
	sub new_gene_entry{	# subroutine to create a new gene entry.
		my $trans = shift;	# reference to transcript hash.
		return {
			'chrom' => $trans->{'chrom'},
			'strand' => $trans->{'strand'},
			'source' => $trans->{'source'},
			'name' => $trans->{'gname'},
			'transcripts' => {
				$trans->{'tid'} => {
					'exons' => $trans->{'exons'}
				}
			},
			'copies' => 1	# number of gene copies.
		};
	}
	sub merge_new_trans{	# subroutine to merge a new transcript into existing gene structure.
		my($gtbl, $trans) = @_;
		my $new_copy_tag = 1;	# Boolean tag for new gene copy. Default=yes.
		my $new_gene_id;	# gene id for the new transcript.
		if($trans->{'chrom'} eq $gtbl->{$trans->{'gid'}}{'chrom'} and 
			$trans->{'strand'} eq $gtbl->{$trans->{'gid'}}{'strand'}){
			$new_copy_tag = 0;
			$new_gene_id = $trans->{'gid'};
		}else{
			for(my $i = 2; $i <= $gtbl->{$trans->{'gid'}}{'copies'}; $i++){
				my $gene_copy_id = $trans->{'gid'} . '#' . $i;
				if($trans->{'chrom'} eq $gtbl->{$gene_copy_id}{'chrom'} and 
					$trans->{'strand'} eq $gtbl->{$gene_copy_id}{'strand'}){
					$new_copy_tag = 0;
					$new_gene_id = $gene_copy_id;
					last;
				}
			}
		}
		if($new_copy_tag){
			# Same gene id, but on different chromosome or strand (a new gene copy).
			$gtbl->{$trans->{'gid'}}{'copies'}++;
			$new_gene_id = $trans->{'gid'} . '#' . $gtbl->{$trans->{'gid'}}{'copies'};
			$gtbl->{$new_gene_id} = new_gene_entry($trans);
		}else{
			# A new transcript within the same gene.
			$gtbl->{$new_gene_id}{transcripts}{$trans->{'tid'}} = {
				'exons' => $trans->{'exons'}
			};
		}
	}
	# Go through each GTF line.
	while(<GTF>){
		chomp;
		my($chrom,$source,$feature,$start,$end,$score,$strand,
			$frame,$attributes) = split /\t/;
		next if $feature ne 'exon';	# ignore non-exon info.
		my %attr_table = ($attributes =~ /(\w+) \"(\S+)\"\;?/g);
		my $gene_id = (exists $attr_table{gene_id})? $attr_table{gene_id} : '';
		my $gene_name = (exists $attr_table{gene_name})? $attr_table{gene_name} : '';
		my $transcript_id = (exists $attr_table{transcript_id})? $attr_table{transcript_id} : '';
		if($chrom eq $cur_transcript{'chrom'} and $strand eq $cur_transcript{'strand'} and
			$gene_id eq $cur_transcript{'gid'} and $gene_name eq $cur_transcript{'gname'} and
			$transcript_id eq $cur_transcript{'tid'}){	# add new exon to transcript structure.
			push @{$cur_transcript{'exons'}},
				{'start' => $start, 'end' => $end, 'class' => 'canonical'};
		}else{
			# Add current transcript to gene table.
			if(exists $gene_table{$cur_transcript{'gid'}}){
				merge_new_trans(\%gene_table, \%cur_transcript);
			}elsif($cur_transcript{'chrom'} ne ''){	# a new non-empty gene id.
				$gene_table{$cur_transcript{'gid'}} = new_gene_entry(\%cur_transcript);
			}
			# Create entry for the new transcript.
			%cur_transcript = (	# working transcript.
				'chrom' => $chrom,
				'strand' => $strand,
				'source' => $source,
				'gid' => $gene_id,
				'gname' => $gene_name,
				'tid' => $transcript_id,
				'exons' => [{'start' => $start, 'end' => $end, 'class' => 'canonical'}],
			);
		}
	}
	# Add last transcript to gene table.
	if(exists $gene_table{$cur_transcript{'gid'}}){
		merge_new_trans(\%gene_table, \%cur_transcript);
	}else{	# a new gene id.
		$gene_table{$cur_transcript{'gid'}} = new_gene_entry(\%cur_transcript);
	}
	close GTF;
	return %gene_table;
}

# Annotate the whole gene table and classify exons.
sub anno_gene_table{
	my $gene_table = shift;	# reference to gene hash table.
	# Annotate each exon and classify them as: promoter, polyA, variant, canonical and altBoundary.
	while(my($gene_id, $gene_struct) = each %{$gene_table}){
		my @trans_ids = keys %{$gene_struct->{transcripts}};
		# Sort exons by start position first.
		foreach my $transcript(@trans_ids){
			my @sorted_exons = sort {$a->{'start'} <=> $b->{'start'}} @{$gene_struct->{transcripts}{$transcript}{exons}};
			$gene_struct->{transcripts}{$transcript}{exons} = [@sorted_exons];
		}
		# Perform pairwise comparison for all transcripts.
		my $strand = $gene_struct->{strand};
		if(@trans_ids > 1){
			for my $i(0..$#trans_ids-1){
				for my $j($i..$#trans_ids){
					my $idA = $trans_ids[$i];
					my $idB = $trans_ids[$j];
					&anno_exon($gene_struct->{transcripts}{$idA}{exons}, $gene_struct->{transcripts}{$idB}{exons}, $strand);
					&anno_exon($gene_struct->{transcripts}{$idB}{exons}, $gene_struct->{transcripts}{$idA}{exons}, $strand);
				}
			}
		}else{
			my $idA = $trans_ids[0];
			&anno_exon($gene_struct->{transcripts}{$idA}{exons}, $gene_struct->{transcripts}{$idA}{exons}, $strand);
		}
	}
}

# Function for exon annotation.
sub anno_exon{
	my($exons_ref,$exons_working,$strand) = @_;
	my $ref_exon_n = @{$exons_ref};
	my $working_exon_n = @{$exons_working};
	my $w = 0;	# counter for working exon.
	my $r = 0;	# iterator for reference exon.
	## Function to determine whether two exons overlap.
	sub overlapExon{
		my($e1,$e2) = @_;
		return ($e1->{start} <= $e2->{end} and $e1->{end} >= $e2->{start});
	}
	foreach my $working_exon(@{$exons_working}){
		$w++;
		# The first and last working exons are directly assigned class.
		if(($w==1 and $strand eq '+') or ($w==$working_exon_n and $strand eq '-')){
			$working_exon->{class} = 'promoter';
		}elsif(($w==1 and $strand eq '-') or ($w==$working_exon_n and $strand eq '+')){
			$working_exon->{class} = 'polyA';
		}else{	# deal with working exon in the middle.
			next if $working_exon->{class} eq 'variant';	# "variant" has priority over other types.
			# Mutually exclusive conditions...
			while($r < $ref_exon_n and $working_exon->{end} >= $exons_ref->[$r]{start}){
				$r++;
			}
			next if $r==0;	# working exon is before/after the reference promoter/polyA.
			next if ($r==$ref_exon_n and $working_exon->{start} > $exons_ref->[$r-1]{end});	# similar as above.
			if(!overlapExon($working_exon, $exons_ref->[$r-1])){	# variant exon.
				$working_exon->{class} = 'variant';
			}else{	# overlap: alternative boundaries.
				next if $working_exon->{class} eq 'altBoth';	# "altBoth" has 2nd priority after "variant".
				my $tentative_class;	# classification based on current overlapping info.
				if($r==1 and $working_exon->{end} != $exons_ref->[$r-1]{end}){	# overlapping left-most ref exon.
					$tentative_class = $strand eq '+'? 'altDonor' : 'altAcceptor';
				}elsif($r==$ref_exon_n and $working_exon->{start} != $exons_ref->[$r-1]{start}){	# overlapping right-most ref exon.
					$tentative_class = $strand eq '+'? 'altAcceptor' : 'altDonor';
				}elsif($r > 1 and $r < $ref_exon_n){	# overlapping middle ref exon.
					if($working_exon->{start} != $exons_ref->[$r-1]{start} and 
						$working_exon->{end} == $exons_ref->[$r-1]{end}){
						$tentative_class = $strand eq '+'? 'altAcceptor' : 'altDonor';
					}elsif($working_exon->{start} == $exons_ref->[$r-1]{start} and 
						$working_exon->{end} != $exons_ref->[$r-1]{end}){
						$tentative_class = $strand eq '+'? 'altDonor' : 'altAcceptor';
					}elsif($working_exon->{start} != $exons_ref->[$r-1]{start} and
						$working_exon->{end} != $exons_ref->[$r-1]{end}){
						$tentative_class = 'altBoth';
					}
				}
				# Must consider previous class and current assignment.
				if(defined $tentative_class){
					if(($working_exon->{class} eq 'altDonor' and $tentative_class eq 'altAcceptor') or 
						($working_exon->{class} eq 'altAcceptor' and $tentative_class eq 'altDonor')){
						$working_exon->{class} = 'altBoth';	# "upgrade" to both alternatives.
					}else{
						$working_exon->{class} = $tentative_class;	# tentative class prevails.
					}
				}
			}
		}
	}
}

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

MyBioinfo::GTF - Read and parse a GTF (Gene Transfer File).

=head1 SYNOPSIS

  use MyBioinfo::GTF;


=head1 DESCRIPTION

  Read a GTF into memory as a hash table. Parse the gene hash table to classify exons.

=head2 EXPORT

  read_gtf_tbl
  anno_gene_table

=head1 SEE ALSO

  Nothing else.

=head1 AUTHOR

Li Shen, E<lt>li.shen@mssm.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2011 by Li Shen

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.6.0 or,
at your option, any later version of Perl 5 you may have available.


=cut
