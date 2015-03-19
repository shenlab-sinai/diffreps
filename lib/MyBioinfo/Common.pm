package MyBioinfo::Common;

use 5.006;
use strict;
use warnings;
use constant INFINITE => 'Inf';
use constant EPSILON => 1e-8;
use POSIX qw(floor ceil);
use Math::CDF;
require MyShortRead::SRBed;
use Data::Dumper;
require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use MyBioinfo::Common ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
	
) ] );

our @EXPORT_OK = qw(mean_r mad padjBH raw_sum2 raw_sum_mean raw_sum_var MchooseN BH_fdr raw_sum_dir raw_sum nb_pval_v2 nb_pval raw_mean_dir raw_mean nb_stat var fold_change chi_stat readnamelist readnamewithinfolist array2hash max min sum mean median log2 log10 read_norm2 rescale_cutoff read_cutoff isAboveCutoff rescale_norm_max rescale_norm_sum1 is_all_zero fprecision unique);

our @EXPORT = qw(padjBH fold_change max min sum mean geomean var median log2 log10 read_norm2 is_all_zero fprecision INFINITE EPSILON);

our $VERSION = '0.61';


######## Preloaded methods go here. ##############

# Imitate the R unique function.
sub unique{
	my %ut;
	foreach(@_){
		$ut{$_} = 1;
	}
	return keys %ut;
}

# Given a vector of P-values, return the adjusted P-values 
# according to the BH procedure.
sub padjBH {
	my $p = shift;	# reference to P-value vector.
	my $n;	# number of tests to multiply.
	if(@_ > 0) {$n = shift;}
	else {$n = @{$p};}
	my %p_BH;	# hash to adjusted P-values.
	# Store indices and P-values into hash.
	for my $i(0..$#{$p}){
		$p_BH{$i} = $p->[$i];
	}
	# Find the sorted indices by raw P-values. This determines the ranks of the P-values.
	my @sorted_r = sort {$p_BH{$a} <=> $p_BH{$b}} keys %p_BH;
	# Apply BH formula to all P-values. 
	# Create a reverse hash from the sorted indices to raw P-value rank at the same time.
	my %rev_r;
	for my $r(1..@sorted_r){
		my $f = $p_BH{$sorted_r[$r-1]} * $n / $r;	# BH formula.
		$p_BH{$sorted_r[$r-1]} = $f > 1? 1.0 : $f;	# truncate values larger than 1.0
		$rev_r{$sorted_r[$r-1]} = $r-1;	# hash: original position -> raw P-value rank.
	}
	# Now, find the sorted indices by BH'ed P-values.
	my @sorted_rr = sort {$p_BH{$a} <=> $p_BH{$b}} keys %p_BH;
	# Go through the 2nd list of sorted indices, solve the inconsistent P-values.
	my $sta_r = 0;	# remember starting raw P-value rank to be adjusted. 
	my $raw_r = 0;	# current raw P-value rank.
	for my $i(0..$#sorted_rr){
		# Sequ: Iterator -> original position -> raw P-value rank.
		if($rev_r{$sorted_rr[$i]} > $raw_r){	# Found a new min P-value. Or else, it is already adjusted.
			$raw_r = $rev_r{$sorted_rr[$i]};	# update current rank.
			for my $j($sta_r..$raw_r-1){
				# use the raw P-value rank to find the original position to be fixed.
				$p_BH{$sorted_r[$j]} = $p_BH{$sorted_rr[$i]};	
			}
			$sta_r = $raw_r + 1;	# advance from current rank.
		}
	}
	# Return the adjusted P-values.
	my @padj;
	for my $i(0..$#{$p}) {push @padj, $p_BH{$i};}
	return @padj;
}

sub raw_sum_mean{
	my ($m1,$m2,$norm_ref)=@_;
	my $v = 0;	
	for my $i(@{$norm_ref})
		{
			$v +=1/$i;
		}		
	
	return $v*($m1+$m2)/2;	
}

sub raw_sum_var{
	my ($base_var,$m1,$m2,$norm,$raw_mean,$eps)=@_;
	my $base_mean = ($m1+$m2)/2;
	my $z =0;
	my $s = 0;
	for my $i(@{$norm})
	{
		$z +=$i;
		$s += (1/$i)**2;
	}
	$z = $z/@{$norm};	
	
	my $var = $base_var - $z * $base_mean;
		if($var < $eps*$base_mean)
		{
			$var = $eps*$base_mean;
		}
		
	return $var*$s + $raw_mean;
	
}


sub MchooseN
{
	my ($m,$n,$ref) = @_;
	my @items = @{$ref};
	
	my $k = $m-$n;
	my @res = ();
	if($n==0)
	{
		return @res;
	}
	if($n==1)
	{
		foreach my $v(@items)
		{
			my @val = ();
			push @val,$v;
			push @res,[@val];
		}
		return @res;
	}
	else
	{
		#to avoid the duplicated combination, treat the items as ordered
		for my $i(0..$k)
		{				
			my @left=();
			for my $j($i+1..$#items)
			{
				push @left,$items[$j];
			}			
			my @ret = MchooseN($m-1,$n-1,\@left);
			my $len = @ret;
		    		    
			for(my $j=0;$j<$len;$j++)
			{
				my @val = ();				
				push @val, $items[$i];
				push @val, @{$ret[$j]};
				
				push @res,[@val];
			}
		}
		return @res;
	}
}

sub BH_fdr{
	my ($p,$c,$threshold,$res) = @_;
	my $index = 1;
	my $max_id =0;
	my $min = 1.0;
	
	foreach my $key (sort {$p->{$a} <=> $p->{$b}} keys %{$p} )
	{
		my $v = $p->{$key};
		$v = $v * $c/$index; 
		
		$res->{$key} = $v;
		
		if($v<=$threshold)
		{
			$res->{$key} = $v;
			if($max_id < $index)
			{
		 		$max_id = $index;
			}
			if($v<$min)
			{
				$min = $v;
			}
		}
		$index ++;
	}
	
	foreach my $key (sort {$p->{$a} cmp $p->{$b}} keys %{$p} )
	{
		if($max_id>0)
		{
			$res->{$key} = $min;
			
			$max_id--;
		}
	}
}

#mean read counts
sub raw_sum{
	my($rarr_srbed,$chrom,$i) = @_;
	my @arr_read;
	# retrieve read count at window# 'i'.
	my $r_mu = 0;
	foreach my $b(@{$rarr_srbed}) {$r_mu += $b->get_bin_count($chrom,$i);}
	
	return $r_mu;
}

sub raw_sum2{
	my($rarr_srbed,$norm_ref) = @_;
	my @arr_read;
	# retrieve read count at window# 'i'.
	my $r_mu = 0;
	foreach my $i(0..@{$rarr_srbed}-1) {$r_mu += $rarr_srbed->[$i]/$norm_ref->[$i];}
	
	return $r_mu;
}
sub raw_sum_dir{
	my($rarr_srbed,$chrom,$i) = @_;
	my @arr_read;
	# retrieve read count at window# 'i'.
	my $r_mu = 0;
	foreach my $b(@{$rarr_srbed}) 
	{
		my @count_both= $b->get_win_count_direction($chrom,$i);
		$r_mu +=$count_both[0]+$count_both[1];
	}
	
	return $r_mu;
}

sub raw_mean{
	my($rarr_srbed,$chrom,$i) = @_;
	my @arr_read;
	# retrieve read count at window# 'i'.
	my $r_mu = 0;
	foreach my $b(@{$rarr_srbed}) {$r_mu += $b->get_bin_count($chrom,$i);}
	
	return $r_mu/@{$rarr_srbed};
}
sub raw_mean2{
	my($m1,$m2,$norm_ref) = @_;
	
	my $v = 0;	
	for my $i(@{$norm_ref})
		{
			$v +=1/$i;
		}	
	$v = $v/@{$norm_ref};
	return $v*($m1+$m2)/2;
}

sub raw_mean_dir{
	my($rarr_srbed,$chrom,$i) = @_;
	my @arr_read;
	# retrieve read count at window# 'i'.
	my $r_mu = 0;
	foreach my $b(@{$rarr_srbed}) {
		my @count_both = $b->get_win_count_direction($chrom,$i);
		$r_mu += $count_both[0] + $count_both[1]; 
		#$b->get_bin_count($chrom,$i);
	}
	
	return $r_mu/@{$rarr_srbed};
}

# Calculate fold change given two values.
sub fold_change{
	my($t,$c) = @_;
	if($t < 0 or $c < 0){
		warn "Negative value in fold change calculation!\n";
		return 0;
	}
	if($t >= $c){
		if($c == 0) {return $t+1;}
		else {return $t/$c;}
	}
	else{
		if($t == 0) {return -($c+1);}
		else {return -$c/$t;}
	}
}

# Calculate Pearson's Chi-square test statistic.
sub chi_stat{
	my($o,$n,$p) = @_;
	return ($o - $n*$p)**2 / ($n*$p*(1-$p));
}

# Function to read a list of names. Assuming the 1st column.
# Default: convert all names to upper case.
# assume the name is in the 1st column and there is no whitespace in names.
sub readnamelist{
  my($file, $ra) = @_;
  my $flag = 1;
  $flag = $_[2] if @_ > 2; # read flag for upper case conversion.
  @{$ra} = (); # deplete the list array first.
  open HLIST, "<", $file or die "Open input file error: $!\n";
  while(<HLIST>){
    chomp;
    my($name) = split;
    $name = uc $name if $flag;
    push @{$ra}, $name;
  }
  close HLIST;
}

# Function to read a list of names with information in 1st and 2nd columns.
# Default: convert all names to upper case.
# assume there is no whitespace in names.
sub readnamewithinfolist{
  my($file, $rh) = @_;
  my $flag = 1;
  $flag = $_[2] if @_ > 2; # read flag for upper case conversion.
  %{$rh} = (); # deplete the list array first.
  open HLIST, "<", $file or die "Open input file error: $!\n";
  while(<HLIST>){
    chomp;
    my($name, $info) = split;
    $name = uc $name if $flag;
    $rh->{$name} = $info;
  }
  close HLIST;
}

# Function to convert an array of names to hash table.
sub array2hash{
  my($ra, $rh) = @_;
  %{$rh} = (); # deplete the hash table first.
  foreach(@{$ra}){
    $rh->{$_} = 1;
  }
}

# Function to find the max element for an array.
sub max{
	die "max function called for an empty array!\n" if @_ < 1;
	my $m = $_[0];
	for(my $i = 1; $i < @_; $i++) {$m = $_[$i] if $_[$i] > $m;}
	return $m;
}

# Function to find the min element for an array.
sub min{
	die "min function called for an empty array!\n" if @_ < 1;
	my $m = $_[0];
	for(my $i = 1; $i < @_; $i++) {$m = $_[$i] if $_[$i] < $m;}
	return $m;
}

# Function to find the sum for an array.
sub sum{
	die "sum function called for an empty array!\n" if @_ < 1;
	my $s = 0;
	foreach my $n(@_) {$s += $n;}
	return $s;
}

# Function to find the mean for an array.
sub mean{
	die "mean function called for an empty array!\n" if @_ < 1;
	my $m = 0;
	my $c = 0;
	foreach my $n(@_) {$m += $n; $c++;}
	return $m / $c;
}

# Function to perform trimmed mean. The array is passed as a reference.
sub mean_r{
	my $r_v = shift;
	if(@{$r_v} == 0){
		warn "Empty array encountered! Return zero.\n";
		return 0;
	}
	my($left_trim, $right_trim);
	if(@_ > 0){	# left trim parameter.
		$left_trim = shift;
		unless($left_trim >= 0 and $left_trim < 0.5){
			warn "Trimming parameter must be in: [0, 0.5). Reset to zero!\n";
			$left_trim = 0;
		}
	}else{
		$left_trim = 0;
	}
	if(@_ > 0){	# right trim parameter.
		$right_trim = shift;
		unless($right_trim >= 0 and $right_trim < 0.5){
			warn "Trimming parameter must be in: [0, 0.5). Reset to zero!\n";
			$right_trim = 0;
		}
	}else{
		$right_trim = 0;
	}
	my @sorted_v = sort {$a <=> $b} @{$r_v};
	my $left_mark = floor($left_trim * @sorted_v);	# start position.
	my $right_mark = @sorted_v - floor($right_trim * @sorted_v);	# end +1 position.
	my $m = 0;
	for(my $i = $left_mark; $i < $right_mark; $i++){
		$m += $sorted_v[$i];
	}
	$m /= $right_mark - $left_mark;
	return $m;
}

# Function to calculate the geometric mean for an array.
sub geomean{
	die "geomean function called for an empty array!\n" if @_ < 1;
	my $m = 0;
	my $c = 0;
	foreach my $n(@_) {
		if($n == 0) {return 0;}
		$m += log($n); 
		$c++;
	}
	return exp($m/$c);
}

# Function to find the variance for an array.
sub var{
	die "variance function called for an empty array!\n" if @_ < 1;
	my $m = 0;
	my $c = 0;
	my $mean = mean(@_);
	foreach my $n(@_) {$m += ($n-$mean)**2; $c++;}
	return $m / ($c-1);
}

# Function to find the median of a numeric array.
sub median{
	die "median function called for an empty array!\n" if @_ < 1;
	my @sorted = sort {$a <=> $b} @_;
	if(@sorted % 2 == 1) {return $sorted[int(@sorted/2)];}
	else {return ($sorted[int(@sorted/2)]+$sorted[int(@sorted/2)-1])/2;}
}

# Function calculate mad: Median absolute deviation.
sub mad{
	my $r_v = shift;
	if(@{$r_v} == 0){	# reference to array.
		return 0;
	}
	my($center, $constant);
	if(@_ > 0){	# center.
		$center = shift;
	}else{
		$center = median(@{$r_v});
	}
	if(@_ > 0){	# scale constant.
		$constant = shift;
		unless($constant > 0){
			warn "Scale constant must be a positive number. Use default!\n";
			$constant = 1.4826;
		}
	}else{
		$constant = 1.4826;
	}
	my @dev;
	foreach my $n(@{$r_v}){
		push @dev, abs($n - $center);
	}
	return $constant * median(@dev);
}

# logorithm base 2.
sub log2{
	my $n = shift;
	if($n == 0) {return -1*INFINITE;}
	return log($n) / log(2);
}

# logorithm base 10.
sub log10{
	my $n = shift;
	if($n == 0) {return -1*INFINITE;}
	return log($n) / log(10);
}

# A subroutine to read normalization constants for treatment and control.
# Syntax:	treatment norm1 norm2...[normN]
# 			control norm1 norm2...[normN]
# whitespace should be used as field separator.
# only identifier 'treatment' and 'control' are recognized and they are case-sensitive.
# only the first two lines of the text file are considered and the rest are ignored.
sub read_norm2{
	my($nf,$rt,$rc) = @_;
	open HNORM, "<", $nf or die "Error in reading $nf:$!\n";
	my @buf = <HNORM>;	# read in all lines into buffer.
	chomp @buf;
	my $tag_t = 0;
	my $tag_c = 0;
	# only deal with the first two lines.
	my @line1 = split ' ', $buf[0];
	if(@line1 > 0 and $line1[0] eq 'treatment'){
		for(my $i = 1; $i < @line1; $i++) {push @{$rt}, $line1[$i];}
		$tag_t = 1;
	}
	elsif(@line1 > 0 and $line1[0] eq 'control'){
		for(my $i = 1; $i < @line1; $i++) {push @{$rc}, $line1[$i];}
		$tag_c = 1;
	}
	my @line2 = split ' ', $buf[1];
	if(@line2 > 0 and $line2[0] eq 'treatment'){
		for(my $i = 1; $i < @line2; $i++) {push @{$rt}, $line2[$i];}
		$tag_t = 1;
	}
	elsif(@line2 > 0 and $line2[0] eq 'control'){
		for(my $i = 1; $i < @line2; $i++) {push @{$rc}, $line2[$i];}
		$tag_c = 1;
	}
	return ($tag_t and $tag_c);	# indicate whether both conditions are met.
}


# A subroutine to read cutoff thresholds for all samples.
# Syntax:	treatment cutoff1 cutoff2...[cutoffN]
# 			control cutoff1 cutoff2...[cutoffN]
# whitespace should be used as field separator.
# only identifier 'treatment' and 'control' are recognized and they are case-sensitive.
# only the first two lines of the text file are considered and the rest are ignored.
sub read_cutoff{
	my($nf,$tr_cut,$co_cut) = @_;
	open HNORM, "<", $nf or die "Error in reading $nf:$!\n";
	my @buf = <HNORM>;	# read in all lines into buffer.
	chomp @buf;
	my $tag_t = 0;
	my $tag_c = 0;
	# only deal with the first two lines.
	my @line1 = split ' ', $buf[0];
	if(@line1 > 0 and $line1[0] eq 'treatment'){
		for(my $i = 1; $i < @line1; $i++) {push @{$tr_cut}, $line1[$i];}
		$tag_t = 1;
	}
	elsif(@line1 > 0 and $line1[0] eq 'control'){
		for(my $i = 1; $i < @line1; $i++) {push @{$co_cut}, $line1[$i];}
		$tag_c = 1;
	}
	my @line2 = split ' ', $buf[1];
	if(@line2 > 0 and $line2[0] eq 'treatment'){
		for(my $i = 1; $i < @line2; $i++) {push @{$tr_cut}, $line2[$i];}
		$tag_t = 1;
	}
	elsif(@line2 > 0 and $line2[0] eq 'control'){
		for(my $i = 1; $i < @line2; $i++) {push @{$co_cut}, $line2[$i];}
		$tag_c = 1;
	}
	return ($tag_t and $tag_c);	# indicate whether both conditions are met.
}

#check if the read counts of current bin is above the cutoff
sub isAboveCutoff
{
	my ($rftr_read,$rfco_read,$rftr_cut,$rfco_cut)=@_;
	my $len = @{$rftr_read};
	
	my $valid = 0;
	
	for my $i(0..@{$rftr_read}-1)
	{
		if($rftr_read->[$i] >= $rftr_cut->[$i])
		{
			$valid = 1;
			
		}
	}	
	if($valid)
	{
		return $valid;
	}	
		for my $i(0..@{$rfco_read}-1)
		{
			if($rfco_read->[$i] >= $rfco_cut->[$i])
			{
				return 1;				
			}
		}
	
	return 0;
}

# A subroutine to rescale normalization constants to the maximal one.
# max_n should be the maximum of both treatment and control.
# However, we may rescale treatment and control separately.
sub rescale_norm_max{
	my($rn,$max_n) = @_;
	my $i = 0;
	foreach my $n(@{$rn}){
		if($n <= 0){
			warn "Normalization constant must be larger than zero! Force it to be zero.\n";
			$rn->[$i++] = 0;
		}
		else {$rn->[$i++] = $max_n / $n;}
	}
}

# A subroutine to rescale cutoff constants.
sub rescale_cutoff{
	my($r_cut,$r_norm) = @_;		
	
	for my $i(0..@{$r_cut}-1){
		$r_cut->[$i] = $r_cut->[$i] * $r_norm->[$i];
	}
}

# A subroutine to rescale normalization constants so that they sum up to one. 
# sum_n should be the summation of both treatment and control.
# However, we may rescale treatment and control separately.
sub rescale_norm_sum1{
	my($rn,$sum_n) = @_;
	my $i = 0;
	foreach my $n(@{$rn}){
		if($n <= 0){
			warn "Normalization constant must be larger than zero! Force it to be zero.\n";
			$rn->[$i++] = 0;
		}
		else {$rn->[$i++] = $n / $sum_n;}
	}
}

# A subroutine to determine whether an array contains all zero elements.
sub is_all_zero{
	my $ra = shift;
	foreach (@{$ra}) {return 0 if $_ != 0;}
	return 1;
}

# Format a scalar or an array of numbers to specified decimal number.
sub fprecision{
	return if @_ < 2;
	my $n = shift;
	if(@_ == 1) {return sprintf "%.$n" . "f", $_[0];}
	else{
		my @a;
		foreach(@_) {push @a, sprintf "%.$n" . "f", $_;}
		return @a;
	}
}
#compute weight factor used in estimating variance
sub compute_beta{
	my ($I,$v,$s,$S) = @_;
	my $beta = (2*($I-1)/($v+2))*(1/$I+$s*$s/$S);
	if($beta>1)
	{return 1;}
	return $beta;
}
#adjust variance
sub adj_var{
	my ($b,$rep_s,$neib_s) = @_;
	
	return (1-$b)*$rep_s+$b*$neib_s;
}
#compute negative binomial statistic score (estimate variance using the variances of its neighbors)
sub nb_stat{
	my ($ref_q,$q_size,$start_pos,$num_rep,$epsilon,$step,$cur_pos) = @_;
		
	#my $cur_pos = $q_size/2;
	if($start_pos ==0 )
	{
		$start_pos = $q_size;
		$cur_pos = 0;
	}
	if($start_pos ==$cur_pos-1)
	{
		$cur_pos =0;
	}
	for(my $l=$start_pos;$l>=$cur_pos;$l--) #process all upstream windows of the window pointed by $cur_pos
	{
		
		my %cur_stat= %{$ref_q->[$l]};
		# do statistics on normalized tr_read and co_read arrays if pass cutoff.
		
		my $num_neighbor = 0;									
		#variance of neighbors
		for(my $k=$q_size;$k>=0;) 
		{			
			my %stat = %{$ref_q->[$k]};	
			$cur_stat{tr_neighbor_var} += $stat{tr_replicate_var}; 
			$cur_stat{co_neighbor_var} += $stat{co_replicate_var};
			$k = $k - $step;
			$num_neighbor++;
		}						
							
		$cur_stat{tr_neighbor_var} = $cur_stat{tr_neighbor_var}/$num_neighbor;	
		$cur_stat{co_neighbor_var} = $cur_stat{co_neighbor_var}/$num_neighbor;							
										
		for(my $k=$q_size;$k>=0;) #variance diff
		{
			my %stat = %{$ref_q->[$k]};
			$cur_stat{tr_var_diff} += ($stat{tr_replicate_var}-$cur_stat{tr_neighbor_var})**2; 	
			$cur_stat{co_var_diff} += ($stat{co_replicate_var}-$cur_stat{co_neighbor_var})**2; 
			$k = $k - $step;							
		}								
							
		my $tr_beta = compute_beta($num_neighbor,$num_rep,$cur_stat{tr_neighbor_var},$cur_stat{tr_var_diff});
		my $co_beta = compute_beta($num_neighbor,$num_rep,$cur_stat{co_neighbor_var},$cur_stat{co_var_diff});	
	#	print "$tr_beta\t$co_beta\n";
		
		my $tr_var = adj_var($tr_beta,$cur_stat{tr_replicate_var},$cur_stat{tr_neighbor_var});
		my $co_var = adj_var($co_beta,$cur_stat{co_replicate_var},$cur_stat{co_neighbor_var});
							
		my $tr_mean = mean(@{$cur_stat{tr_read}});
		my $co_mean = mean(@{$cur_stat{co_read}});
							
		if($tr_mean > $co_mean) {$cur_stat{dirn}='Up';}
		if($tr_mean < $co_mean) {$cur_stat{dirn}='Down';}
		if($tr_mean ==$co_mean) {$cur_stat{dirn}='--';}				
							
		if($tr_var <$tr_mean)
		{
			$tr_var = $tr_mean + $epsilon;
		}
		if($co_var < $co_mean)
		{
			$co_var = $co_mean + $epsilon;
		}
				
	#	print "$cur_stat{tr_mean}\t$tr_var\t$cur_stat{co_mean}\t$co_var\n";			
	    $tr_var = $cur_stat{tr_replicate_var};
	    $co_var = $cur_stat{co_replicate_var};
		my $pval = nb_pval($cur_stat{tr_mean},$cur_stat{co_mean},$tr_mean,$tr_var,$co_mean,$co_var,$epsilon);
		
		$cur_stat{score} = $pval;
	#	print "$cur_stat{tr_mean}\t$cur_stat{tr_replicate_var}\t$tr_var\t$cur_stat{co_mean}\t$cur_stat{co_replicate_var}\t$co_var\n";											
	#	print "pval = $pval\n";
		
		$ref_q->[$l] = \%cur_stat;
							
	}						
}

sub nb_pval_v2{
	my ($ka,$miu1,$var1,$eps) = @_;
	
	
	my @rp1 = ();
	nb_r_p($miu1,$var1,\@rp1,$eps);
		
	my $r1 = $rp1[0];
	my $p1 = $rp1[1];
		
	if($ka <= $miu1){
		return &Math::CDF::pnbinom($ka,$r1,$p1);
	}
	else
	{
		return 1 - &Math::CDF::pnbinom($ka,$r1,$p1);
	}		
}

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

MyBioinfo::Common - Commonly used Perl subroutines for my bioinformatics work.

=head1 SYNOPSIS

  use MyBioinfo::Common;
  or
  use MyBioinfo::Common qw(mean_r mad padjBH raw_sum2 raw_sum_mean raw_sum_var MchooseN BH_fdr raw_sum_dir raw_sum nb_pval_v2 nb_pval raw_mean_dir raw_mean nb_stat var fold_change chi_stat readnamelist readnamewithinfolist array2hash max min sum mean median log2 log10 read_norm2 rescale_cutoff read_cutoff isAboveCutoff rescale_norm_max rescale_norm_sum1 is_all_zero fprecision unique)

  my $s = &nb_stat($ref_q, $q_size, $start_pos, $num_rep, $epsilon, $step, $cur_pos);
  my $fc = &fold_change($treatment_val, $control_val);
  my $s = &chi_stat($obs, $n, $prob);
  &readnamelist($filename, \@array);
  &readnamewithinfolist($filename, \@hash);
  &array2hash(\@ref_array, \@ref_hash);
  my $m = &max(@array);
  my $m = &min(@array);
  my $s = &sum(@array);
  my $m = &mean(@array);
  my $m = &median(@array);
  my $v = &log2($m);
  &read_norm2($filename, \@treatment, \@control);
  &rescale_norm_max(\@normalization, $max_norm);
  &rescale_norm_sum1(\@normalization, $sum_norm);
  my $boolean = &is_all_zero(\@array);
  my @formatted = &fprecision($n_after_decimal, @array);
  &read_cutoff($filename, \@treatment, \@control);

=head1 DESCRIPTION

  Some convenient functions for my bioinformatics work. 


=head2 EXPORT

  fold_change
  chi_stat
  readnamelist
  readnamewithinfolist
  array2hash
  max
  min
  sum
  mean
  median
  log2
  read_norm2
  rescale_norm_max
  rescale_norm_sum1
  is_all_zero
  fprecision


=head1 SEE ALSO

MyBioinfo::Math

MyBioinfo::NBTest

=head1 AUTHOR

Li Shen, E<lt>shenli.sam@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2010-2013 by Li Shen

diffReps goes under GNU GPL v3: http://www.gnu.org/licenses/gpl.html


=cut
