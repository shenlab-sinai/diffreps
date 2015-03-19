use 5.006;
our $VERSION = '1.122';

# Class for a chromatin modification site consisting consecutive significant windows.
package diffReps::ChromaModSite;
use strict;
use MyBioinfo::Common;
use MyShortRead::SRBed;

our @ISA = qw(diffReps::SlideWindow);	# Inherited from 'SlideWindow'.

# Constructor.
sub new{
	my $class = shift;
	my $self = {
		chrom => undef,
		start => undef,
		end => undef,
		siteCenter => undef,
		retrieved => 0,		# boolean tag: whether count data have been retrieved?
		tr_cnt => [],		# region normalized treatment counts.
		co_cnt => [],		# region normalized control counts.
		tr_cnt_ => [],		# normalized treatment counts using #reads.
		co_cnt_ => [],		# normalized control counts using #reads.
		btr_cnt => [],		# normalized background treatment counts.
		bco_cnt => [],		# normalized background control counts.
		#trEnr => undef,		# treatment enrichment vs. background.
		#coEnr => undef,		# control enrichment vs. background.
		rsumTr => undef,	# sum of raw treatment counts.
		rsumCo => undef,	# sum of raw control counts.
		dirn => undef,		# change direction: 'Up' or 'Down'.
		logFC => undef,		# log2 fold change.
		pval => undef,		# P-value.
		padj => undef,		# adjusted P-value.
		winSta => undef,	# best window start.
		#winTr => [],		# window normalized treatment counts.
		#winCo => [],		# window normalized control counts.
		winFC => undef,		# window log2 fold change.
		winP => undef,		# window P-value.
		winQ => undef		# window adjusted P-value.
	};
	bless $self, $class;
	return $self;
}

# Copy member data into an anonymous hash table.
sub dcopy{
	my $self = shift;
	my $new = {};
	%{$new} = %{$self};
	# Deep copy array members that are array variables.
	$new->{tr_cnt} = []; @{$new->{tr_cnt}} = @{$self->{tr_cnt}};
	$new->{co_cnt} = []; @{$new->{co_cnt}} = @{$self->{co_cnt}};
	$new->{tr_cnt_} = []; @{$new->{tr_cnt_}} = @{$self->{tr_cnt_}};
	$new->{co_cnt_} = []; @{$new->{co_cnt_}} = @{$self->{co_cnt_}};
	$new->{btr_cnt} = []; @{$new->{btr_cnt}} = @{$self->{btr_cnt}};
	$new->{bco_cnt} = []; @{$new->{bco_cnt}} = @{$self->{bco_cnt}};
	#$new->{winTr} = []; @{$new->{winTr}} = @{$self->{winTr}};
	#$new->{winCo} = []; @{$new->{winCo}} = @{$self->{winCo}};
	bless $new, 'diffReps::ChromaModSite';
	return $new;
}

# Initialize a region by a significant window.
sub init{
	my($self,$gdt,$rwin) = @_;
	$self->{chrom} = $rwin->{chrom};
	$self->{start} = $rwin->{winN}*$gdt->{step} + 1;
	$self->{end} = $rwin->{winN}*$gdt->{step} + $gdt->{winSize};
	$self->{siteCenter} = ($self->{start} + $self->{end}) /2;
	$self->{retrieved} = 1;
	@{$self->{tr_cnt}} = @{$rwin->{tr_cnt}};
	@{$self->{co_cnt}} = @{$rwin->{co_cnt}};
	@{$self->{tr_cnt_}} = @{$rwin->{tr_cnt_}};
	@{$self->{co_cnt_}} = @{$rwin->{co_cnt_}};
	@{$self->{btr_cnt}} = @{$rwin->{btr_cnt}};
	@{$self->{bco_cnt}} = @{$rwin->{bco_cnt}};
	$self->{rsumTr} = $rwin->{rsumTr};
	$self->{rsumCo} = $rwin->{rsumCo};
	$self->{dirn} = $rwin->{dirn};
	$self->{logFC} = $rwin->{logFC};
	$self->{pval} = $rwin->{pval};
	$self->{winSta} = $self->{start};
	#@{$self->{winTr}} = @{$rwin->{tr_cnt}};
	#@{$self->{winCo}} = @{$rwin->{co_cnt}};
	$self->{winFC} = $rwin->{logFC};
	$self->{winP} = $rwin->{pval};
}

# Expand the genomic coordinates. Replace best window if necessary.
# Do not retrieve counts and sums to save some computation.
sub expand{
	my($self,$gdt,$rwin) = @_;
	$self->{end} = $rwin->{winN}*$gdt->{step} + $gdt->{winSize};	# new region end.
	$self->{siteCenter} = ($self->{start} + $self->{end}) /2;
	# Invalidate previous diff results.
	$self->reset_content(1);	# reset everything but keep direction.
#	$self->{tr_cnt} = [];
#	$self->{co_cnt} = [];
#	$self->{rsumTr} = undef;
#	$self->{rsumCo} = undef;
#	$self->{logFC} = undef;
	# Replace best window if necessary.
    if($rwin->{pval} < $self->{winP}){
    	$self->{winSta} = $rwin->{winN}*$gdt->{step} + 1;
    	#@{$self->{winTr}} = @{$rwin->{tr_cnt}};
    	#@{$self->{winCo}} = @{$rwin->{co_cnt}};
    	$self->{winFC} = $rwin->{logFC};
    	$self->{winP} = $rwin->{pval};
    }
}

# Print a formatted record of the current region. Do nothing if the P-val is not defined.
sub print_reg{
	my($self,$gdt,$h) = @_;
    if(defined $self->{pval}){
    	# Concatenate treatment and control counts for output.
    	my @tr_cnt_fmt = fprecision(2, @{$self->{tr_cnt}});
    	my @co_cnt_fmt = fprecision(2, @{$self->{co_cnt}});
    	my $tr_cnt_str = join(';', @tr_cnt_fmt);
    	my $co_cnt_str = join(';', @co_cnt_fmt);
		my($tr_bkg_enr, $co_bkg_enr);
		if($gdt->{bkgEnr}){
			my $btr_m = mean(@{$self->{btr_cnt}});
			my $bco_m = mean(@{$self->{bco_cnt}});
			$tr_bkg_enr = $btr_m > 0? fprecision(2, mean(@{$self->{tr_cnt_}}) / $btr_m) : INFINITE;
			$co_bkg_enr = $bco_m > 0? fprecision(2, mean(@{$self->{co_cnt_}}) / $bco_m) : INFINITE;
		}else{
			$tr_bkg_enr = 'NA';
			$co_bkg_enr = 'NA';
		}
    	# Output: chrom,start,end,length,(treatment read count),(control read count),direction,logFC,pval,padj,
    	# (window start),(window end),(window logFC),(window pval),(window padj).
    	print $h join("\t", ($self->{chrom},
        	$self->{start}, $self->{end}, $self->{end}-$self->{start}+1,	# region coordinates.
        	$tr_cnt_str, $co_cnt_str,	# formatted read counts.
			fprecision(2, &mean(@{$self->{tr_cnt}})), fprecision(2, &mean(@{$self->{co_cnt}})),	# avg read count.
			$tr_bkg_enr, $co_bkg_enr,	# enrichment vs. background.
        	$self->{dirn}, fprecision(2,$self->{logFC}), # change direction and logFC.
    		$self->{pval}, $self->{padj},	# diff analysis info.
        	$self->{winSta}, $self->{winSta}+$gdt->{winSize}-1,	# best window coordinates.
        	fprecision(2,$self->{winFC}), 
    		$self->{winP}, $self->{winQ})), "\n";	# best window diff info.
    }
}

# Do we need to perform stat test for the region?
# If the region contains only one window, this is already done.
sub needStat{
	my $self = shift;
	return (defined($self->{winP}) && !defined($self->{pval}));
}

# Retrieve normalized read counts for the region. Override parent function.
sub retrieve_norm_cnt{
	my($self,$gdt) = @_;
	my @atr_cnt;	# array treatment counts.
	my $i = 0;	# iterator for normalization constants.
	foreach my $b(@{$gdt->{bedTr}}) {
		push @atr_cnt, ($b->get_region_count($self->{chrom}, $self->{start}, $self->{end}) 
			/ $gdt->{normTr}[$i++]);
	}
	my @aco_cnt;	# array control counts.
	$i = 0;
	foreach my $b(@{$gdt->{bedCo}}) {
		push @aco_cnt, ($b->get_region_count($self->{chrom}, $self->{start}, $self->{end}) 
			/ $gdt->{normCo}[$i++]);
	}
	@{$self->{tr_cnt}} = @atr_cnt;
	@{$self->{co_cnt}} = @aco_cnt;
	# Retrieve norm count for background if needed.
	if($gdt->{bkgEnr}){
		my @atr_cnt_;	# array treatment counts.
		$i = 0;	# iterator for normalization constants.
		foreach my $b(@{$gdt->{bedTr}}) {
			push @atr_cnt_, ($b->get_region_count($self->{chrom}, $self->{start}, $self->{end}) / $gdt->{normTr_}->[$i++]);
		}
		@{$self->{tr_cnt_}} = @atr_cnt_;
		my @aco_cnt_;	# array treatment counts.
		$i = 0;	# iterator for normalization constants.
		foreach my $b(@{$gdt->{bedCo}}) {
			push @aco_cnt_, ($b->get_region_count($self->{chrom}, $self->{start}, $self->{end}) / $gdt->{normCo_}->[$i++]);
		}
		@{$self->{co_cnt_}} = @aco_cnt_;
		if(@{$gdt->{bedBtr}}){
			my @abtr_cnt;	# array background treatment counts.
			$i = 0;	# iterator for normalization constants.
			foreach my $b(@{$gdt->{bedBtr}}) {
				push @abtr_cnt, ($b->get_region_count($self->{chrom}, $self->{start}, $self->{end}) / $gdt->{normBtr}->[$i++]);
			}
			@{$self->{btr_cnt}} = @abtr_cnt;
		}
		if(@{$gdt->{bedBco}}){
			my @abco_cnt;	# array background treatment counts.
			$i = 0;	# iterator for normalization constants.
			foreach my $b(@{$gdt->{bedBco}}) {
				push @abco_cnt, ($b->get_region_count($self->{chrom}, $self->{start}, $self->{end}) / $gdt->{normBco}->[$i++]);
			}
			@{$self->{bco_cnt}} = @abco_cnt;
		}
		if(@{$self->{btr_cnt}} == 0){
			@{$self->{btr_cnt}} = @{$self->{bco_cnt}};
		}
		if(@{$self->{bco_cnt}} == 0){
			@{$self->{bco_cnt}} = @{$self->{btr_cnt}};
		}
	}
}

# Retrieve sum of raw read counts for the region. Override parent function.
sub retrieve_raw_sum{
	my($self,$gdt) = @_;
	$self->{rsumTr} = 0;
	foreach my $b(@{$gdt->{bedTr}}) {
		$self->{rsumTr} += $b->get_region_count($self->{chrom}, $self->{start}, $self->{end});
	}
	$self->{rsumCo} = 0;
	foreach my $b(@{$gdt->{bedCo}}) {
		$self->{rsumCo} += $b->get_region_count($self->{chrom}, $self->{start}, $self->{end});
	}
}

# Determine whether a window is within gap.
sub isWithinGap{
	my($self,$gdt,$rwin) = @_;
	return 0 if !defined $self->{chrom};
	my $winSta = $rwin->{winN}*$gdt->{step}+1;
	return ($winSta - $self->{end}-1 <= $gdt->{gap} and $rwin->{dirn} eq $self->{dirn});
}


# Class for a sliding window of fixed size.
package diffReps::SlideWindow;

use strict;
use Statistics::TTest;
use MyBioinfo::NBTest;
use MyBioinfo::Common;
use MyBioinfo::Math qw( gtest_gof chisqtest_gof );
use MyShortRead::SRBed;

# Constructor.
sub new{
	my $class = shift;
	my $self = {
		chrom => undef,
		winN => undef,		# window number to retrieve count from SRBed class.
		maxN => undef,		# maximum window number.	
		retrieved => 0,		# boolean tag: whether count data have been retrieved?
		tr_cnt => [],		# normalized treatment counts.
		co_cnt => [],		# normalized control counts.
		tr_cnt_ => [],		# normalized treatment counts using #reads.
		co_cnt_ => [],		# normalized control counts using #reads.
		btr_cnt => [],		# normalized background treatment counts.
		bco_cnt => [],		# normalized background control counts.
		rsumTr => undef,	# sum of raw treatment counts.
		rsumCo => undef,	# sum of raw control counts.
		dirn => undef,		# change direction: 'Up' or 'Down'.
		logFC => undef,		# log2 fold change.
		pval => undef		# P-value.
	};
	bless $self, $class;
	return $self;
}

# Set the current chromosome name, max window number and reset everything.
# Obtain counts and sums for the first window.
sub set_chrom{
	my($self,$gdt,$chrom,$maxN) = @_;
	$self->{chrom} = $chrom;
	$self->{maxN} = $maxN;
	$self->{winN} = 0;
	$self->reset_content;
	#$self->{retrieved} = 0;
	#	$self->retrieve_norm_cnt($gdt);
	#	$self->retrieve_raw_sum($gdt);
	#$self->{trEnr} = undef;
	#$self->{coEnr} = undef;
#	$self->{dirn} = undef;		# change direction: 'Up' or 'Down'.
#	$self->{logFC} = undef;		# log2 fold change.
#	$self->{pval} = undef;		# P-value.
}

# Reset window number.
sub reset_winN{
	my $self = shift;
	$self->{winN} = 0;
}

# Reset window/region content.
sub reset_content{
	my($self,$keepdi) = @_;
	$self->{retrieved} = 0;	# boolean tag: whether count data have been retrieved?
	$self->{tr_cnt} = [];		# normalized treatment counts.
	$self->{co_cnt} = [];		# normalized control counts.
	$self->{tr_cnt_} = [];	# normalized treatment counts using #reads.
	$self->{co_cnt_} = [];	# normalized control counts using #reads.
	$self->{btr_cnt} = [];	# normalized background treatment counts.
	$self->{bco_cnt} = [];	# normalized background control counts.
	$self->{rsumTr} = undef;	# sum of raw treatment counts.
	$self->{rsumCo} = undef;	# sum of raw control counts.
	$self->{dirn} = undef unless $keepdi;	# change direction: 'Up' or 'Down'.
	$self->{logFC} = undef;	# log2 fold change.
	$self->{pval} = undef;		# P-value.
}

# Retrieve norm/raw count given current window/region location.
sub retrieve_new_cnt{
	my($self,$gdt) = @_;
	$self->{retrieved} = 1;
	$self->retrieve_norm_cnt($gdt);
	$self->retrieve_raw_sum($gdt);
}

# Shift the window one position to the right and obtain new counts and sums.
sub move_forward{
	my($self,$gdt) = @_;
	if($self->{winN} <= $self->{maxN}-1) {
		$self->{winN}++;
		$self->reset_content;
		#	$self->retrieve_norm_cnt($gdt);
		#	$self->retrieve_raw_sum($gdt);
		#$self->{trEnr} = undef;
		#$self->{coEnr} = undef;
    	#$self->{dirn} = undef;		# change direction: 'Up' or 'Down'.
    	#$self->{logFC} = undef;		# log2 fold change.
    	#$self->{pval} = undef;		# P-value.
		return $self->{winN};
	}
#	elsif($self->{winN} == $self->{maxN}-1){
#		$self->{winN} = $self->{maxN};
#		$self->{tr_cnt} = [];
#		$self->{co_cnt} = [];
#		$self->{tr_cnt_} = [];
#		$self->{co_cnt_} = [];
#		$self->{btr_cnt} = [];
#		$self->{bco_cnt} = [];
#		#$self->{trEnr} = undef;
#		#$self->{coEnr} = undef;
#    	$self->{dirn} = undef;		# change direction: 'Up' or 'Down'.
#    	$self->{logFC} = undef;		# log2 fold change.
#    	$self->{pval} = undef;		# P-value.
#		return $self->{winN};
#	}
	else {return $self->{maxN};}
	# Max valid window index should be maxN-1.
}

# Return current window number.
sub get_winN{
	my $self = shift;
	if($self->{winN} < $self->{maxN}) {return $self->{winN};}
	else {return -1;}
}

# Return chromosome name.
sub get_chrom{
	my $self = shift;
	return $self->{chrom};
}

# Get normalized read counts for external procedures.
sub get_cnt_tr{
	my $self = shift;
	return @{$self->{tr_cnt}};
}
sub get_cnt_co{
	my $self = shift;
	return @{$self->{co_cnt}};
}
sub get_cnt_tr_{
	my $self = shift;
	return @{$self->{tr_cnt_}};
}
sub get_cnt_co_{
	my $self = shift;
	return @{$self->{co_cnt_}};
}
sub get_cnt_btr{
	my $self = shift;
	return @{$self->{btr_cnt}};
}
sub get_cnt_bco{
	my $self = shift;
	return @{$self->{bco_cnt}};
}

# Get differential analysis info.
sub get_diff_info{
	my $self = shift;
	return {
		dirn =>		$self->{dirn},
		logFC =>	$self->{logFC},
		pval =>		$self->{pval}
	};	# anonymous hash.
}

# Retrieve normalized read counts by current window number from SRBed objects.
sub retrieve_norm_cnt{
	my($self,$gdt) = @_;
	my @atr_cnt;	# array treatment counts.
	my $i = 0;	# iterator for normalization constants.
	foreach my $b(@{$gdt->{bedTr}}) {
		push @atr_cnt, ($b->get_win_count($self->{chrom}, $self->{winN}) / $gdt->{normTr}->[$i++]);
	}
	@{$self->{tr_cnt}} = @atr_cnt;
	my @aco_cnt;	# array control counts.
	$i = 0;
	foreach my $b(@{$gdt->{bedCo}}) {
		push @aco_cnt, ($b->get_win_count($self->{chrom}, $self->{winN}) / $gdt->{normCo}->[$i++]);
	}
	@{$self->{co_cnt}} = @aco_cnt;
	# Retrieve norm count for background if needed.
	if($gdt->{bkgEnr}){
		my @atr_cnt_;	# array treatment counts.
		$i = 0;	# iterator for normalization constants.
		foreach my $b(@{$gdt->{bedTr}}) {
			push @atr_cnt_, ($b->get_win_count($self->{chrom}, $self->{winN}) / $gdt->{normTr_}->[$i++]);
		}
		@{$self->{tr_cnt_}} = @atr_cnt_;
		my @aco_cnt_;	# array treatment counts.
		$i = 0;	# iterator for normalization constants.
		foreach my $b(@{$gdt->{bedCo}}) {
			push @aco_cnt_, ($b->get_win_count($self->{chrom}, $self->{winN}) / $gdt->{normCo_}->[$i++]);
		}
		@{$self->{co_cnt_}} = @aco_cnt_;
		if(@{$gdt->{bedBtr}}){
			my @abtr_cnt;	# array background treatment counts.
			$i = 0;	# iterator for normalization constants.
			foreach my $b(@{$gdt->{bedBtr}}) {
				push @abtr_cnt, ($b->get_win_count($self->{chrom}, $self->{winN}) / $gdt->{normBtr}->[$i++]);
			}
			@{$self->{btr_cnt}} = @abtr_cnt;
		}
		if(@{$gdt->{bedBco}}){
			my @abco_cnt;	# array background treatment counts.
			$i = 0;	# iterator for normalization constants.
			foreach my $b(@{$gdt->{bedBco}}) {
				push @abco_cnt, ($b->get_win_count($self->{chrom}, $self->{winN}) / $gdt->{normBco}->[$i++]);
			}
			@{$self->{bco_cnt}} = @abco_cnt;
		}
		if(@{$self->{btr_cnt}} == 0){
			@{$self->{btr_cnt}} = @{$self->{bco_cnt}};
		}
		if(@{$self->{bco_cnt}} == 0){
			@{$self->{bco_cnt}} = @{$self->{btr_cnt}};
		}
	}
}

# Retrieve sum of raw read counts by current window number from SRBed objects.
sub retrieve_raw_sum{
	my($self,$gdt) = @_;
	$self->{rsumTr} = 0;
	foreach my $b(@{$gdt->{bedTr}}) {
		$self->{rsumTr} += $b->get_win_count($self->{chrom}, $self->{winN});
	}
	$self->{rsumCo} = 0;
	foreach my $b(@{$gdt->{bedCo}}) {
		$self->{rsumCo} += $b->get_win_count($self->{chrom}, $self->{winN});
	}
}

# Determine whether current bin counts are above cutoff or not.
sub isAboveCutoff{
	my($self,$gdt) = @_;
	die "Empty SRBed vectors encountered in function: isAboveCutoff.\n" 
		if @{$gdt->{bedTr}} == 0 or @{$gdt->{bedCo}} == 0;
	my $tagTr = 1;
	my $tagCo = 1;
	if($gdt->{bkgEnr} and $gdt->{useBkg} > 0){	# use background as filter.
		if(!$self->{retrieved}){	# retrieve new count if not done so.
			$self->retrieve_new_cnt($gdt);
		}
		my $btr_m = mean(@{$self->{btr_cnt}});
		my $bco_m = mean(@{$self->{bco_cnt}});
		foreach my $c(@{$self->{tr_cnt_}}){
			if(!$c or ($c and $btr_m and $c / $btr_m <= $gdt->{useBkg})){
				$tagTr = 0;
				last;
			}
		}
		if(!$tagTr){
			foreach my $c(@{$self->{co_cnt_}}){
				if(!$c or ($c and $bco_m and $c / $bco_m <= $gdt->{useBkg})){
					$tagCo = 0;
					last;
				}
			}
		}
	}else{	# use raw count and equation as filter.
		my $i = 0;
		foreach my $b(@{$gdt->{bedTr}}){
			my $rawcnt = $b->get_win_count($self->{chrom}, $self->{winN});
			if(!$rawcnt or $rawcnt <= $gdt->{meanTr}[$i] + $gdt->{nsd}*$gdt->{stdTr}[$i]){
				$tagTr = 0;
				last;
			}
			$i++;
		}
		if(!$tagTr){
			$i = 0;
			foreach my $b(@{$gdt->{bedCo}}){
				my $rawcnt = $b->get_win_count($self->{chrom}, $self->{winN});
				if(!$rawcnt or $rawcnt <= $gdt->{meanCo}[$i] + $gdt->{nsd}*$gdt->{stdCo}[$i]){
					$tagCo = 0;
					last;
				}
				$i++;
			}
		}
	}
	return ($tagTr || $tagCo);
}

# Calculate statistical scores based on current window counts.
sub calc_stat{
	my($self,$gdt,$meth,$eps) = @_;
	#return $self->{pval} if defined $self->{pval};	# no duplicated calculation please.
	return 1 if !defined($self->{chrom});	# do nothing if the region is not defined.
	if(!$self->{retrieved}){	# retrieve new count if not done so.
		$self->retrieve_new_cnt($gdt);
	}
	# Check whether counts and sums are ready for stat test.
	die "Empty count vectors encountered when doing stat test.\n" 
		if @{$self->{tr_cnt}} == 0 or @{$self->{co_cnt}} == 0;
	die "Undefined raw sums for treatment or control when doing stat test.\n"
		if !defined($self->{rsumTr}) or !defined($self->{rsumCo});
	# Noramlized mean count for each group.
	my $tr_m = mean(@{$self->{tr_cnt}});
	my $co_m = mean(@{$self->{co_cnt}});
	# Special case: zero count for both groups.
	if($tr_m == 0 and $co_m == 0) {
		$self->{logFC} = 0;
		$self->{dirn} = 'Nons';
		$self->{pval} = 1;
	}
	else {
		if($co_m == 0) {$self->{logFC} = INFINITE;}
		else {$self->{logFC} = log2($tr_m/$co_m);}
		$self->{dirn} = $self->{logFC} > 0? 'Up' : 'Down';
		# Perform statistical tests to find P-value.
		if($meth eq 'tt'){		# T-test.
			my $ttest = new Statistics::TTest;
			$ttest->load_data(\@{$self->{tr_cnt}}, \@{$self->{co_cnt}});
			$self->{pval} = $ttest->t_prob;
		}
		elsif($meth eq 'nb'){	# Negative Binomial test.
			my $tr_var = var(@{$self->{tr_cnt}});
			my $co_var = var(@{$self->{co_cnt}});

			my $tr_raw_m = raw_sum_mean($tr_m,$co_m,$gdt->{normTr});
			my $co_raw_m = raw_sum_mean($tr_m,$co_m,$gdt->{normCo});
			my $tr_raw_v = raw_sum_var($tr_var,$tr_m,$co_m,$gdt->{normTr},$tr_raw_m,$eps);
			my $co_raw_v = raw_sum_var($co_var,$tr_m,$co_m,$gdt->{normCo},$co_raw_m,$eps);

			$self->{pval} = nb_pval($self->{rsumTr},$self->{rsumCo},$tr_raw_m,$tr_raw_v,$co_raw_m,$co_raw_v,$eps);
		}
		elsif($meth eq 'gt' or $meth eq 'cs'){
			my $tr_sum = sum(@{$self->{tr_cnt}});
			my $co_sum = sum(@{$self->{co_cnt}});
			my $a_cnt = [$tr_sum, $co_sum];	# ref to array of normalized counts.
			my $tr_rep = scalar @{$self->{tr_cnt}};
			my $co_rep = scalar @{$self->{co_cnt}};
			my $tot_rep = $tr_rep + $co_rep;
			my $a_p = [$tr_rep/$tot_rep, $co_rep/$tot_rep];	# ref to array of probabilities.
			my @test_res;	# result vector: statistic, DOF, p-value.
			if($meth eq 'gt'){
				@test_res = gtest_gof($a_cnt, $a_p);
			}else{
				@test_res = chisqtest_gof($a_cnt, $a_p);
			}
			$self->{pval} = $test_res[2];
		}
		else{
			die "Undefined statistical test! Halt execution.\n";
		}
	}
	return $self->{pval};
}


# Subroutines for NB test. Copied from Common.pm. Written by Ying Jin.
# NOTES: In the original diffPermute program, norm constants were multiplied to
# raw counts to get normalized values. This only saves ignorable execution time
# but causes confusion. They are now reversed, so do the following subroutines.
sub raw_sum_mean{
	my ($m1,$m2,$norm_ref)=@_;
	my $v = 0;	
	for my $i(@{$norm_ref}){
			$v += $i;
	}
	
	return $v*($m1+$m2)/2;	
}
sub raw_sum_var{
	my ($base_var,$m1,$m2,$norm,$raw_mean,$eps)=@_;
	my $base_mean = ($m1+$m2)/2;
	my $z =0;
	my $s = 0;
	for my $i(@{$norm}){
		$z += 1/$i;
		$s += $i**2;
	}
	$z = $z/@{$norm};	
	
	my $var = $base_var - $z * $base_mean;
	if($var < $eps*$base_mean){
		$var = $eps*$base_mean;
	}
		
	return $var*$s + $raw_mean;
	
}
###################################



# Class for managing a list of significant regions: add, adjust P-values, output.
package diffReps::RegionList;
use strict;
use MyBioinfo::Common;
use diffReps::DiffRes qw( cmpsite );

# Constructor.
sub new{
	my $class = shift;
	my $self = {
		regList => [],	# region list.
		adjusted => 0	# bool tag for P-value adjustment.
	};
	bless $self, $class;
	return $self;
}

# Add a significant region into list if the P-value is defined.
sub add{
	my($self,$rreg) = @_;
	if(defined $rreg->{pval}){
		push @{$self->{regList}}, $rreg->dcopy();
	}
}

# Add a significant region into list if the P-value is defined.
# Delete the original region.
sub add_del{
	my($self,$rreg) = @_;
	if(defined $rreg->{pval}){
		push @{$self->{regList}}, $rreg->dcopy();
		%{$rreg} = %{new diffReps::ChromaModSite};
	}
}

# Append another region list to the end.
sub append_list{
	my($self,$rl) = @_;
	if(@{$rl->{'regList'}}){
		push @{$self->{'regList'}}, @{$rl->{'regList'}};
		$self->{'adjusted'} = 0;
	}
}

# Sort all differential sites according to genomic coordinates.
sub sort{
	my $self = shift;
	@{$self->{regList}} = sort cmpsite @{$self->{regList}};
}

# Adjust P-values for all regions in the list.
sub adjPval{
	my $self = shift;
	return if $self->{adjusted};
	my $N;	# the 'N' in BH formula.
	if(@_ > 0) {$N = shift;}
	else {$N = @{$self->{regList}};}
	# Extract all P-values from the list and feed them to an adjustment procedure.
	my(@p1, @p2);	# pair of region and window P-values.
	foreach my $r(@{$self->{regList}}) {
		push @p1, $r->{pval};
		push @p2, $r->{winP};
	}
	my @q1 = padjBH(\@p1, $N);
	my @q2 = padjBH(\@p2, $N);
	my $i = 0;	# iterator for adjusted P-value vector.
	foreach my $r(@{$self->{regList}}) {
		$r->{padj} = $q1[$i];
		$r->{winQ} = $q2[$i++];
	}
	# Set bool tag.
	$self->{adjusted} = 1;
}

# Print header line to file handle.
sub gen_header{
	my($self,$hrep) = @_;
	# Print header line.
	print $hrep "Chrom\tStart\tEnd\tLength\tTreatment.cnt\tControl.cnt\tTreatment.avg\tControl.avg\tTreatment.enr\tControl.enr\tEvent\tlog2FC\tpval\tpadj\twinSta\twinEnd\twinFC\twinP\twinQ\n";
}

# Output all formatted regions.
sub output{
	my($self,$gdt,$h) = @_;
	if(!$self->{adjusted}) {$self->adjPval();}
	foreach my $r(@{$self->{regList}}) {$r->print_reg($gdt,$h);}
}




# Preloaded methods go here.

# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

diffReps::SlideWindow - class to deal with a sliding window.

diffReps::ChromaModSite - class to deal with chromatin modification site, inherited from SlideWindow.

diffReps::RegionList - class for a list of modification regions.

=head1 SYNOPSIS

  my $slideWindow = new diffReps::SlideWindow;
  my $site = new diffReps::ChromaModSite;
  my $regList = new diffReps::RegionList;

=head1 DESCRIPTION

In order to perform differential analysis for ChIP-seq data, we need to define classes and methods
that can retrieve short read counts for a designated window; perform statistics and book-keep important
information.

=head2 EXPORT

None. This is OO designed.



=head1 SEE ALSO

diffReps::DiffRes

Mailing list: https://groups.google.com/forum/#!forum/diffreps-discuss

Web site: https://code.google.com/p/diffreps/


=head1 AUTHOR

Li Shen, E<lt>shenli.sam@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010-2013 by Li Shen

diffReps goes under GNU GPL v3: http://www.gnu.org/licenses/gpl.html


=cut
