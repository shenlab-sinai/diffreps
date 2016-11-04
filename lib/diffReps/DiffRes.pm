use 5.006;
our $VERSION = '1.11.2';

# Package to analyze the ChIP-seq differential analysis results.
# Contains functions to manipulate and score a cluster of diff sites.
# Can be used to identify chromatin modifiation hotspots.

package diffReps::DiffRes;
use strict;
#use MyShortRead::MyShortRead qw(compare2);
use Math::CDF qw(ppois);
use MyBioinfo::Common qw( min );
use constant INF => 1e10;
# use constant WIN1 => 2e5;	# 200kb.
# use constant WIN2 => 1e6;	# 1Mb.

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw( cmpsite );


# Constructor.
# Lambda: defined as the avg site-to-site distance which can be used to 
# calculate an expected segment number for a cluster.
# Expected segments = cluster distance / lambda
# P-value is evaluated as Poisson(Observed segment |Expected segment).
# With the same Observed / Expected ratio, larger cluster is always favored 
# vs. smaller one.
sub new{
	my $class = shift;
	my($w1, $w2) = @_;
	my $self = {
		siteList => [],			# diff site list.
		# each site is a hash: (chrom, start, end, siteCenter, recLoc).
		curPos => 0,			# current position in site list.
		glbLambda => undef,		# global lambda: avg inter-distance between two sites.
		win1 => $w1,			# 1st local window size.
		win2 => $w2,			# 2nd local window size.
		clust => {
			from => undef,		# start position (in the list) of cluster.
			to => undef,		# end position (in the list) of cluster.
			seg => undef,		# cluster segments.
			dist => undef,		# left-to-right distance.
			enrich => undef,	# observed / expected segment.
			pval => undef		# P-value based on Poisson distribution.
		}
	};
	bless $self, $class;
	return $self;
}

# Evaluate the significance of current cluster.
sub sigeval{
	my($self) = @_;
	$self->dist_seg_clust();	# calc actual distance, segments for cluster.
	# Trivial case: the new pushed site is on another chromosome or the cluster contains a single site.
	if($self->{'clust'}{'dist'} == INF or $self->{'clust'}{'dist'} == 0){
		$self->{'clust'}{'seg'} = undef;
		$self->{'clust'}{'enrich'} = undef;
		$self->{'clust'}{'pval'} = 1.0;
		return 1.0;
	}
	my($l_d1,$l_s1,$l_d2,$l_s2) = $self->calc_locL(-1);	# search left.
	my($r_d1,$r_s1,$r_d2,$r_s2) = $self->calc_locL(+1);	# search right.
	my $d1 = $l_d1 + $r_d1;
	my $s1 = $l_s1 + $r_s1;
	my $d2 = $l_d2 + $r_d2;
	my $s2 = $l_s2 + $r_s2;

	# Exception: the local lambda estimate may become zero if two sites' start
	# positions coincide with each other. To prevent dividing by zero, set the
	# lambda estimate to 1bp.
	my $L1 = $s1 > 0? $d1 / $s1 : INF;	# local lambda 1.
	$L1 = $L1 == 0? 1 : $L1;
	my $L2 = $s2 > 0? $d2 / $s2 : INF;	# local lambda 2.
	$L2 = $L2 == 0? 1 : $L2;

	my $L_min = &min($L1, $L2, $self->{'glbLambda'});
	if($L_min != INF and defined $self->{'clust'}{'dist'}){
		my $exp_seg = $self->{'clust'}{'dist'} / $L_min;
		$self->{'clust'}{'enrich'} = $self->{'clust'}{'seg'} / $exp_seg;
		$self->{'clust'}{'pval'} = 1 - ppois($self->{'clust'}{'seg'}, $exp_seg);
	}else{
		$self->{'clust'}{'enrich'} = undef;
		$self->{'clust'}{'pval'} = 1.0;
	}

	return $self->{'clust'}{'pval'};
}

# Search either direction of a differential list to calculate local lambda. 
sub calc_locL{
	my($self,$di) = @_;	# di must be: -1 or 1 for left or right.
	# Calc site distance considering direction.
	sub site_dist_di{
		my($di,$w_site,$c_site) = @_;
		return $di < 0? &site_dist($w_site, $c_site) : &site_dist($c_site, $w_site);
	}
	# Calc segment number considering direction.
	sub seg_di{
		my($di,$w_loc,$c_loc) = @_;
		return $di < 0? $c_loc - $w_loc : $w_loc - $c_loc;
	}
	# Move back to the last site within window size and calculate distance and segments 
	# considering direction.
	sub last_dist_seg_di{
		my($self,$di,$w_loc,$c_loc,$c_site) = @_;
		my $w_loc_ = $w_loc - $di;	# recover to the last site that is still on the same chromosome.
		my $w_site = $self->{'siteList'}->[$w_loc_];
		my $d = &site_dist_di($di, $w_site, $c_site);
		my $s = &seg_di($di, $w_loc_, $c_loc);

		return ($d,$s);
	}

	# Setup initial parameter values.
	my $c_loc = $di < 0? $self->{'clust'}{'from'} : $self->{'clust'}{'to'};	# cluster edge location.
	my $c_site = $self->{'siteList'}->[$c_loc];	# cluster edge site.
	my($dist1,$seg1,$dist2,$seg2);	# local distance and #segment 1 and 2.
	$dist1 = $dist2 = $seg1 = $seg2 = 0;
	my $tag1 = 0;	# tag for window1.
	my $nsite = @{$self->{'siteList'}};	# number of total sites.
	my $w_loc;
	for($w_loc = $c_loc + $di; $w_loc >= 0 and $w_loc < $nsite; $w_loc += $di){
		my $w_site = $self->{'siteList'}->[$w_loc];	# window site.
		my $d = &site_dist_di($di, $w_site, $c_site);	# distance between two sites.
		if(!$tag1 and $d > $self->{'win1'}){	# pass window1 for the first time?
			($dist1,$seg1) = $self->last_dist_seg_di($di, $w_loc, $c_loc, $c_site);
			$tag1 = 1;	# set tag for window1.
		}
		if($d > $self->{'win2'}){	# pass window2.
			($dist2,$seg2) = $self->last_dist_seg_di($di, $w_loc, $c_loc, $c_site);
			last;
		}
	}
	if($w_loc < 0 or $w_loc >= $nsite){	# exceeded the begining or end of the diff list.
		my($d,$s) = $self->last_dist_seg_di($di, $w_loc, $c_loc, $c_site);
		if(!$tag1){
			($dist1,$seg1) = ($d, $s);
			$tag1 = 1;
		}
		($dist2,$seg2) = ($d, $s);
	}

	return ($dist1, $seg1, $dist2, $seg2);
}

# Read genomic coordinates of chromatin modification sites from a list
# of files. Record the line number for each diff site.
sub read_modsite{
	my($self,$rfile,$info) = @_;
	# Figure out which columns are used in output.
	# Format eg: 1,2,3-6,7-8,9
	my @icol = ();	# columns.
	if($info ne '0'){
		my @iseg = split /,/, $info;	# segments.
		foreach my $s(@iseg){
			my($sta,$end) = split /-/, $s;
			$end = $sta unless defined $end;
			for my $c($sta..$end){
				push @icol, $c;
			}
		}
	}
	# Start extracting sites from inputs.
	$self->{siteList} = [];	# deplete site list first.
	foreach my $f(@{$rfile}){
		open DIFF, '<', $f or die "Open diff site file error: $!\n";
		my $ln = 0;
		my $site;
		while(<DIFF>){	# Get rid of header lines.
			$ln++;
			last if ($site = &parse_site($_, $f, $ln, \@icol));
		}
		push @{$self->{siteList}}, $site if $site;
		while(<DIFF>){
			$ln++;
			if($site = &parse_site($_, $f, $ln, \@icol)){
				push @{$self->{siteList}}, $site;
			}
			else {die "Un-recognized diff site format in file $f: $ln\n";}
		}
		close DIFF;
	}
	$self->{curPos} = 0;	# reset current position.
}


# Parse a line of diff site file and store it in hash table.
sub parse_site{
	my($rec,$file,$ln,$rcol) = @_;
	my $mark;
	# Try to figure out the histone mark name from the input file...
	if($file =~ /((h[0-9]+k[0-9]+((me[0-9]+)|ac))|(polII)|(pol2))/i) {$mark = $1;}
	else {$mark = $file;}	# if not matching pattern, use file name instead.
	if($rec =~ /^(\S+)\t([0-9]+)\t([0-9]+)/){
		my @cells = split /\t/, $rec;
		# Extract output info from specified columns.
		my @vinfo;
		foreach my $c(@{$rcol}){
			if($c > 0 and $c <= @cells){
				push @vinfo, $cells[$c-1];
			}
		}
		my $sinfo = '';
		if(@vinfo > 0){
			$sinfo = join(',', @vinfo);
		}
		# Return record as a hash table.
		return {
			chrom => $cells[0],  # chromosome name.
			start => $cells[1],
			end => $cells[2],
			siteCenter => ($cells[1] + $cells[2]) / 2,  # diff site center.
			recLoc => $file . ':(' . $sinfo . '):' . $ln,  # record location in file:(info):ln format.
			markType => $mark
		};
	}
	else {return 0;}
}

# Sort diff sites according to genomic coordinates.
sub sort_list{
	my($self) = @_;
	if(@{$self->{siteList}} > 1){
		@{$self->{siteList}} = sort cmpsite @{$self->{siteList}};
	}
}

# Comparing two sites according to their centers.
sub cmpsite{
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

  return $a->{siteCenter} <=> $b->{siteCenter};
}

# Calculate global lambda. Return negative value if it cannot be determined.
sub calc_glbL{
	my($self) = @_;
	if(@{$self->{siteList}} < 2){
		$self->{'glbLambda'} = INF;
		return;
	}
	my $tot_dist = 0;	# total internal distance.
	my $tot_seg = 0;	# total distance segments.
	my $left = 0;
	# Go through each chromosome and do:
	# identify the left and right-most site, calculate the distance and
	# record the number of internal segments; add them to global variables.
	# Finally, calcualte an averaged distance.
	my $l_site = $self->{siteList}->[$left];
	my $r_site;
	my $right;
	for($right = 1; $right < @{$self->{siteList}}; $right++){
		$r_site = $self->{siteList}->[$right];
		if(!same_chr($l_site, $r_site)){
			$r_site = $self->{siteList}->[$right-1];	# the last site on same chromosome.
			my $d = site_dist($l_site, $r_site);
			if($d < INF){
				$tot_dist += $d;
				$tot_seg += $right - $left -1;
			}
			$left = $right;	# move into new chromosome.
			$l_site = $self->{siteList}->[$left];
		}
	}	# Upon exit: left and right sites are always on the same chromosome.
	# Add the last chromosome.
	$tot_dist += site_dist($l_site, $r_site);
	$tot_seg += $right - $left -1;
	$self->{'glbLambda'} = $tot_seg > 0? $tot_dist / $tot_seg : INF;
}

# Determine if two sites are on the same chromosome.
sub same_chr{
	my($l_site,$r_site) = @_;
	return $l_site->{'chrom'} eq $r_site->{'chrom'};
}

# Calcualte distance between two sites. Always return a positive value.
# If Left > Right, stop the program and issue a fatal error.
sub site_dist{
	my($l_site,$r_site) = @_;
	if(!same_chr($l_site, $r_site)){
		return INF;
	}
	if($l_site->{'siteCenter'} > $r_site->{'siteCenter'}){
		die "Left site has larger coordinate than right site when calculating distance. Check program logic. Stop now.\n";
	}
	return $r_site->{'siteCenter'} - $l_site->{'siteCenter'};
}

# Initialize a cluster with the current site.
sub init_clust{
	my($self) = @_;
	$self->{clust}{from} = $self->{curPos};
	$self->{clust}{to} = $self->{curPos};
	$self->{clust}{seg} = undef;
	$self->{clust}{dist} = undef;
	$self->{clust}{enrich} = undef;
	$self->{clust}{pval} = 1.0;

	return $self->{clust}{pval};
}

# Are we at the end of the site list?
sub is_list_end{
	my($self) = @_;
	return $self->{curPos} >= @{$self->{siteList}} - 1;
}

# Shift the leftmost site out of the cluster.
sub shift_clust{
	my($self) = @_;
	if($self->{clust}{from} < $self->{clust}{to}){
		$self->{clust}{from}++;
		return $self->sigeval();
	}else{
		warn "Only one site left in cluster. Shift was not performed.\n";
		return $self->{clust}{pval};
	}
}

# Un-shift the leftmost site into the cluster.
sub unshift_clust{
	my($self) = @_;
	if($self->{clust}{from} > 0){
		$self->{clust}{from}--;
		return $self->sigeval();
	}else{
		warn "Leftmost site is already at the beginning of the difflist. Unshift was not performed.\n";
		return $self->{clust}{pval};
	}
}

# Add a site to the rightmost of the cluster.
sub push_clust{
	my($self) = @_;
	if($self->{curPos} < @{$self->{siteList}}){
		$self->{curPos}++;	# move current position to the right.
		$self->{clust}{to}++;
		return $self->sigeval();
	}
	else {
		warn "Current position is already at the end of the list. Push was not performed.\n";
		return $self->{clust}{pval}; 
	}
}

# Remove the rightmost site from the cluster.
sub pop_clust{
	my($self) = @_;
	if($self->{clust}{to} > $self->{clust}{from}){
		$self->{clust}{to}--;
		return $self->sigeval();
	}
	else{
		warn "Site list cannot be empty. Pop was not performed.\n";
		return $self->{clust}{pval};
	}
}

# Update distance and segments for current cluster.
sub dist_seg_clust{
	my($self) = @_;
	my $l_loc = $self->{clust}{from};
	my $r_loc = $self->{clust}{to};
	$self->{'clust'}{'seg'} = $r_loc - $l_loc;
	my $l_site = $self->{'siteList'}->[$l_loc];
	my $r_site = $self->{'siteList'}->[$r_loc];
	$self->{'clust'}{'dist'} = &site_dist($l_site, $r_site);
}

# Extract a cluster's info to external procedure.
sub get_clustinfo{
	my($self,$h) = @_;
	my $from = $self->{clust}{from};
	my $to = $self->{clust}{to};
	if($from < $to){
		my @v_rec;	# vector of records of cluster members.
		for my $i($from..$to){
			push @v_rec, $self->{siteList}[$i]{recLoc};
		}
		my %h_mark;	# hash of mark types;
		for my $i($from..$to){
			$h_mark{$self->{siteList}[$i]{markType}} = 1;
		}
		my @v_mark = sort {$a cmp $b} keys %h_mark;	# vector of unique mark types.
		my $chrom = $self->{siteList}[$from]{chrom};
		my $start = $self->{siteList}[$from]{start};
		my $end = $self->{siteList}[$to]{end};
		return {	# return an anonymous hash.
			chrom => $chrom,
			start => $start,
			end => $end,
			len => $end-$start+1,
			enrich => $self->{clust}{enrich},
			pval => $self->{clust}{pval},
			v_rec => [@v_rec],
			v_mark => [@v_mark]
		};
	}else{
		warn "Cluster contains less than two sites. No print was done.\n";
		return {};
	}
}

# Get number of sites.
sub get_siteN{
	my($self) = @_;
	return scalar @{$self->{siteList}};
}


# Class to manage cluster list.
package diffReps::ClustList;
use strict;
use MyBioinfo::Common qw(padjBH);

# Constructor.
sub new{
	my $class = shift;
	my $self = {
		clustList => [],
		f_handler => undef
	};
	bless $self, $class;
	return $self;
}

# Add a cluster to the list.
sub add{
	my($self,$rc) = @_;
	push @{$self->{clustList}}, $rc;
}

# Adjust cluster P-values by BH.
sub adjPval{
	my($self,$N) = @_;
	my @p;
	foreach my $c(@{$self->{clustList}}) {push @p, $c->{pval};}
	my @padj = padjBH(\@p, $N);
	my $i = 0;	# iterator for adjusted P-values.
	foreach my $c(@{$self->{clustList}}) {$c->{padj} = $padj[$i++];}
}

# Output all clusters' info to a file.
sub output{
	my($self) = @_;
	if(!defined $self->{f_handler}){
		warn "File handler has not been opened yet.\n";
		return;
	}
	my $h = $self->{f_handler};
	foreach my $c(@{$self->{clustList}}){
		print $h join("\t", $c->{chrom}, 
			$c->{start}, $c->{end}, $c->{len},
			$c->{enrich}, 
			$c->{pval}, $c->{padj}, 
			scalar @{$c->{v_rec}}, join(';', @{$c->{v_rec}}),
			scalar @{$c->{v_mark}}, join(';', @{$c->{v_mark}})), "\n";
	}
}

# Output header line.
sub print_header{
	my($self) = @_;
	my $h = $self->{f_handler};
	if(defined $h){
		print $h "Chrom\tStart\tEnd\tLength\tenrich\tpval\tpadj\tnsite\tSites\tntype\tMarkType\n";
	}else{
		warn "File handler has not been opened yet. Nothing was done.\n";
	}
}

# Spit a line into opened file.
sub spit{
	my($self,$line) = @_;
	my $h = $self->{f_handler};
	if(defined $h){
		print $h "$line\n";
	}else{
		warn "File handler has not been opened yet. Nothing was done.\n";
	}
}

# Open file handler.
sub open_file{
	my($self,$f) = @_;
	my $h;
	open $h, '>', $f or die "Open output file error: $!\n";
	$self->{f_handler} = $h;
}

# Close file handler.
sub done{
	my($self) = @_;
	close $self->{f_handler} if defined $self->{f_handler};
}




# Preloaded methods go here.

# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

diffReps::DiffRes - class to store and manipulate differential sites.

diffReps::ClustList - class to deal with a cluster of differential sites(hotspots).

=head1 SYNOPSIS

  my $diffRes = new diffReps::DiffRes;
  my $clustList = new diffReps::ClustList;

=head1 DESCRIPTION

Classes and methods to read differential sites from text files and manipulate them.

=head2 EXPORT

None. This is OO designed.



=head1 SEE ALSO

diffReps::ChromaModSite

Mailing list: https://groups.google.com/forum/#!forum/diffreps-discuss

Web site: https://code.google.com/p/diffreps/


=head1 AUTHOR

Li Shen, E<lt>shenli.sam@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010-2013 by Li Shen

diffReps goes under GNU GPL v3: http://www.gnu.org/licenses/gpl.html


=cut
