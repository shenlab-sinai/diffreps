#!/usr/bin/perl -w
#
# Perform sliding window differential analysis to identify chromatin
# modification sites for ChIP-seq data from two groups of samples.
# Assume the two groups are treatment and control.
#
# Adapted from previous slide_window.pl and re-designed in OOP fashion.
#
# by Li Shen, Dec 2010.
#
# TODO:	1. Add directionality index? (Do we still need this?)
#
# Added: G-test for two groups without replicates. (04/19/2011)
# Changes: re-implemented Chisq test and G-test using my own codes.
#

use 5.006;
use strict;
our $VERSION = '1.55.6';
use Getopt::Long;
use MyBioinfo::Common;
use MyBioinfo::Common qw(mean_r mad);
use MyShortRead::MyShortRead;
use MyShortRead::SRBed;
use diffReps::ChromaModSite;
use PJ::Genome qw( get_chrsz );
use POSIX qw(ceil floor);
use constant GENOM_SAMN => 100000;
use IO::Handle;
use Parallel::ForkManager 0.7.6;
use constant DEBUG => 0;


my @treatment;	# array of treatment file names.
my @control;	# array of control file names.
my @btr;	# background treatment samples.
my @bco;	# background control samples.
my $gname = '';	# genome name.
my($nsd,$sign,$meth,$step,$window,$gap,$chrlen,$norm,$nrpass,$report,$fragsize);
my $std = 0;	# bool tag for using standard estimation.
my $alpha = 0.05;	# alpha for right-trimmed mean.
my $noanno = 0;	# tag for switching off annotation.
my $nohs = 0;	# tag for switching off looking for hotspots.
$nsd = 2;		# default number of std above mean as cutoff.
my $bkg = 0;	# tag for using background as filter. 
$sign = 1e-4;		# default P-value cutoff.
$meth = 'nb';		# default stat test: NB.
$norm = undef;		# norm constants file name.
$nrpass = 1;		# switch: norm on bins passed std cutoff.
$step = -1;		# default step/bin size.
$window = -1;
$gap = 0;		# maximal gap allowed between two continuous regions.
$fragsize = 100;	# default fragment size.
my $nproc = 1;	# number of processors to use.
my $mode = 'peak';	# scanning mode: peak, nucleosome, block.
my $fullcmd = join(' ', $0, @ARGV);	# record command arguments.

if(@ARGV < 1 or !GetOptions('treatment=s{1,}' => \@treatment,
							'control=s{1,}' => \@control,
							'btr:s{1,}' => \@btr,
							'bco:s{1,}' => \@bco,
							'gname:s' => \$gname,
							'chrlen:s' => \$chrlen,
							'report=s' => \$report,
							'mode:s' => \$mode,
							'meth:s' => \$meth,
							'pval:f' => \$sign,
							'norm:s' => \$norm,
							'nrpass:i' => \$nrpass,
							'nsd:s' => \$nsd,
							'bkg:f' => \$bkg,
							'std' => \$std,
							'noanno' => \$noanno,
							'nohs' => \$nohs,
							'alpha:f' => \$alpha,
							'window:i' => \$window,
							'gap:i' => \$gap,
							'step:i' => \$step,
							'nproc:i' => \$nproc,
							'fragsize:i' => \$fragsize)){
	print "\ndiffReps - Detect Differential Sites from ChIP-seq with Biological Replicates.\n\n";
	print "Usage: diffReps.pl --treatment bed_file1...[bed_fileN] --control bed_file1...[bed_fileN] --report output_results \\\n";
	print "           [--gname genome_name|--chrlen chrom_length_file]\n\n";
	print "**Chromosome lengths can be specified through a text file or given a genome name.\n";
	print "**Currently built-in genomes: mm9, hg19, rn4.\n\n";
	print "Optional parameters(defaults in parentheses):\n";
	print "- Background samples(DNA input or IgG control):\n";
	print "    --btr            Treatment group background: bed_file1...[bed_fileN].\n";
	print "    --bco            Control group background: bed_file1...[bed_fileN].\n";
	print "    Hint: If background is only specified for one group, it will automatically be used for both groups.\n\n";
	print "- Genomic region parameters:\n";
	print "    --mode(peak)     Scanning mode: a selection implies a different window size.\n";
	print "                     Set window and step size manually to override.\n";
	print "                     (p)eak      (=1000)  Histone mark peak (Default).\n";
	print "                     (n)ucleosome(=200)   Single nucleosome (+DNAlinker).\n";
	print "                     (b)lock     (=10000) Large chromatin modification block.\n";
	print "    --window(1000)   Window size (default=Histone mark peak size).\n";
	print "    --step(1/10 win) Window moving step size.\n";
	print "    --gap(0)         Gap allowed between two consecutive windows.\n\n";
	print "- Background filtering using: mean + nsd*deviation.\n";
	print "    --std            Use standard estimation of mean and deviation (Default=Robust estimation).\n";
	print "                     In robust estimation, median absolute deviation is used in place of standard deviation.\n";
	print "    --nsd(broad)     Z-score cutoff for low read count. Choose from two default modes or set your own.\n";
	print "                     (b)road     (=2)   Broad peak such as H3K36me3.\n";
	print "                     (s)harp     (=20)  Sharp peak such as H3K4me3 or Transcription factors.\n";
	print "    --alpha(0.05)    Alpha for right-trimmed mean, must be in: [0, 0.5).\n";
	print "    --bkg(0)         Use fold enrichment vs. background as filter instead. Set a float number such as 2.0 here.\n";
	print "                     Default is to use the Z-score as filter.\n\n";
	print "- Statistical testing parameters:\n";
	print "    --meth(nb)       Statistical test (nb=Negative binomial; gt=G-test; tt=T-test; cs=Chi-square test).\n";
	print "    --pval(0.0001)   P-value cutoff for significant windows.\n\n";
	print "- Normalization can be done externally and be supplied as a text file:\n";
	print "    --nrpass(1)      Do normalization on bins pass nsd cutoff?\n";
	print "    --norm           File name to specify pre-determined norm constants (Default=Estimate by diffReps).\n\n";
	print "- Misc. parameters:\n";
	print "    --frag(100)      ChIP-seq library fragment size. Use to shift read positions.\n";
	print "    --nproc(1)       Number of processors to use.\n";
	print "    --noanno         Switch off genomic annotation for differential sites (Default=Do annotation).\n";
	print "    --nohs           Switch off looking for chromatin modification hotspots (Default=Find hotspots).\n";
	exit 0;
}

# Choice of statistical test.
die "Statistic not defined. Exit.\n" unless $meth eq 'tt' or $meth eq 'nb' or $meth eq 'gt' or $meth eq 'cs';
my $statDesc;
if($meth eq 'tt') {$statDesc = 'T-test';}
elsif($meth eq 'nb') {$statDesc = "Negative Binomial";}
elsif($meth eq 'gt') {$statDesc = "G-test";}
elsif($meth eq 'cs') {$statDesc = "Chisquare test";}

if(($meth eq 'tt' or $meth eq 'nb') and (@treatment < 2 or @control < 2)) {
	print "To use ".$statDesc.", you must have at least two replicates per condition.\n";
	print "Use G-test(preferred) or Chisquare test instead. Exit.\n";
	exit 0;
}

# Choice of window and step size.
($window, $step) = choose_winstep($mode, $window, $step);

# Background filtering parameters;
my $robust = {'std' => $std, 'alpha' => $alpha};

# Choose Z-score cutoff.
if($nsd eq 'broad' or $nsd eq 'b'){
	$nsd = 2;
}elsif($nsd eq 'sharp' or $nsd eq 's'){
	$nsd = 20;
}elsif(!($nsd =~ /^[0-9]+$/)){
	warn "Unspecified mode for --nsd. Use default value of 2.\n";
	$nsd = 2;
}

# Package-wide forkmanager.
my $pm = new Parallel::ForkManager($nproc); 

# Open report file for output.
my $hrep;	# report file handle.
open $hrep, ">", $report or die "Cannot create report file:$!\n";

# Record some parameter settings.
spit($hrep, "# diffReps version", $VERSION);
spit($hrep, '## Parameter Settings ##');
spit($hrep, '# Treatment files', @treatment);
spit($hrep, '# Control files', @control);
spit($hrep, '# Treatment background', @btr? @btr : 'None Specified');
spit($hrep, '# Control background', @bco? @bco : 'None Specified');
spit($hrep, '# Output file', $report);
spit($hrep, '# Statistical Test', $statDesc);
spit($hrep, '# Use robust estimation', $std ? 'no' : 'yes');
spit($hrep, '# Number of deviation cutoff', $nsd);
spit($hrep, '# Alpha in right-trimmed mean', $alpha);
spit($hrep, '# P-value cutoff', $sign);
spit($hrep, '# Window size', $window);
spit($hrep, '# Step size', $step);
spit($hrep, '# Gap size', $gap);
spit($hrep, '# Fragment size', $fragsize);
spit($hrep, '# Number of processors used', $nproc);
spit($hrep, '# Full command used', $fullcmd);
$hrep->flush;

# Step0: read chromosome length.
my @chrlen_ordered;
my %chrlen_tbl;
my $ref_chrtbl = get_chrsz($gname);
if(keys %{$ref_chrtbl} == 0){
	if(defined $chrlen){
		@chrlen_ordered = read_chrlen_ordered($chrlen);	# ordered array.
		%chrlen_tbl = read_chrlen_tbl($chrlen);	# hash table.
		$gname = '';	# reset for tagging purpose.
	}else{
		die "Neither correct genome name nor chromosome length file provided. Abort.\n";
	}
}else{
	%chrlen_tbl = %{$ref_chrtbl};
	@chrlen_ordered = get_ordered_chrlen($ref_chrtbl);
}

# Step1: create treatment and control SRBed objects;
# separate BED file into trunks of chromosomes.
my $time_tag = time;
my @tr_srbed;
my @co_srbed;
my @btr_srbed;
my @bco_srbed;
my @all_nread;  # array for all read counts.

if(DEBUG) {  # for debug, do not use forkmanager.
	foreach my $bedfile(@treatment) {
		my $srbed = new MyShortRead::SRBed($bedfile,$window,$fragsize);
		my $nr = $srbed->create_chr_list(\%chrlen_tbl);
		push @tr_srbed, $srbed;
		push @all_nread, $nr;		
	}
	foreach my $bedfile(@control) {
		my $srbed = new MyShortRead::SRBed($bedfile,$window,$fragsize);
		my $nr = $srbed->create_chr_list(\%chrlen_tbl);
		push @co_srbed, $srbed;
		push @all_nread, $nr;		
	}
	if(@btr) {
		foreach my $bedfile(@btr) {
			my $srbed = new MyShortRead::SRBed($bedfile,$window,$fragsize);
			my $nr = $srbed->create_chr_list(\%chrlen_tbl);
			push @btr_srbed, $srbed;
			push @all_nread, $nr;		
		}
	}
	if(@bco) {
		foreach my $bedfile(@bco) {
			my $srbed = new MyShortRead::SRBed($bedfile,$window,$fragsize);
			my $nr = $srbed->create_chr_list(\%chrlen_tbl);
			push @bco_srbed, $srbed;
			push @all_nread, $nr;		
		}
	}
} else {  # production code.
	my %all_nread;  # hash table for all read counts.
	my %all_srbed;	# hash table for all SRBed obj.
	# Create forkmanager for file splitting.
	$pm->run_on_finish(	# callback function to retrieve SRBed objects.
		sub{
			my($pid, $exit_code, $ident, $exit_signal, $core_dump, $dr) = @_;
			if(defined $dr){
				$all_nread{$dr->{'bed'}} = $dr->{'nread'};
				$all_srbed{$dr->{'bed'}} = $dr->{'srbed'};
			}else{
				die "Child process returned void data in splitting BED files with exit code: $exit_code.\n";
			}
		}
	);
	# Generate forks for splitting.
	foreach my $bedfile( (@treatment, @control, @btr, @bco) ){
		my $pid = $pm->start and next;
		my $srbed = new MyShortRead::SRBed($bedfile,$window,$fragsize);
		my $nr = $srbed->create_chr_list(\%chrlen_tbl);
		my %rd = ('bed' => $bedfile, 'nread' => $nr, 'srbed' => $srbed);
		$pm->finish(0, \%rd);
	}
	$pm->wait_all_children;

	# Re-assemble data in hash table into arrays.
	foreach my $b(@treatment){
		push @tr_srbed, $all_srbed{$b};
		push @all_nread, $all_nread{$b};
	}
	foreach my $b(@control){
		push @co_srbed, $all_srbed{$b};
		push @all_nread, $all_nread{$b};
	}
	if(@btr){
		foreach my $b(@btr){
			push @btr_srbed, $all_srbed{$b};
			push @all_nread, $all_nread{$b};
		}
	}
	if(@bco){
		foreach my $b(@bco){
			push @bco_srbed, $all_srbed{$b};
			push @all_nread, $all_nread{$b};
		}
	}
}
spit($hrep, '# Split BED files done in ' . (time - $time_tag) . ' seconds.');
$hrep->flush;

# Step2: find mean, std and normalization constants for all samples.
$time_tag = time;
my @all_beds = (@tr_srbed, @co_srbed);
my $msn_res;	# result of mean, std and norm.
if($nrpass){
	$msn_res = mean_std_norm($robust, \@all_beds, \@chrlen_ordered, $nsd);
}
else {
	$msn_res = mean_std_norm($robust, \@all_beds, \@chrlen_ordered);
}
spit($hrep, '# Estimate mean and deviation for read counts done in ' . (time - $time_tag) . ' seconds');
my(@tr_mean, @co_mean, @tr_std, @co_std);
split_array($msn_res->{mean}, scalar @tr_srbed, \@tr_mean, \@co_mean);
split_array($msn_res->{std}, scalar @tr_srbed, \@tr_std, \@co_std);
spit($hrep, '## Important values estimated from data ##');
spit($hrep, '# Read count statistics based on bin size', $window);
spit($hrep, '# Treatment count mean', fprecision(2, @tr_mean));
spit($hrep, '# Control count mean', fprecision(2, @co_mean));
spit($hrep, '# Treatment count deviation', fprecision(2, @tr_std));
spit($hrep, '# Control count deviation', fprecision(2, @co_std));
my @tr_norm;	# treatment normalization constants.
my @co_norm;	# control normalization constants.
if(defined $norm){
	die "Read normalization constants failed. Check syntax.\n" 
		unless read_norm2($norm, \@tr_norm, \@co_norm);
	die "Number of normalization constants does not equal number of samples.\n" 
		unless @tr_norm == @treatment and @co_norm == @control;
}
else {split_array($msn_res->{norm}, scalar @tr_srbed, \@tr_norm, \@co_norm);}
spit($hrep, '# Treatment normalization constant', fprecision(2, @tr_norm));
spit($hrep, '# Control normalization constant', fprecision(2, @co_norm));
$hrep->flush;
# Do normalization using total read number including background samples.
my @all_norm;
my $all_gm = geomean(@all_nread);
foreach my $nr(@all_nread){
	push @all_norm, $nr / $all_gm;
}
my @tr_norm_ = splice @all_norm, 0, @treatment;
my @co_norm_ = splice @all_norm, 0, @control;
my @btr_norm = splice @all_norm, 0, @btr if @btr;
my @bco_norm = splice @all_norm, 0, @bco if @bco;
# Reset bin and window sizes for sliding window procedure.
foreach my $b(@tr_srbed){
	$b->set_bin_size($step);
	$b->set_window_size($window);
}
foreach my $b(@co_srbed){
	$b->set_bin_size($step);
	$b->set_window_size($window);
}
if(@btr_srbed){
	foreach my $b(@btr_srbed){
		$b->set_bin_size($step);
		$b->set_window_size($window);
	}
}
if(@bco_srbed){
	foreach my $b(@bco_srbed){
		$b->set_bin_size($step);
		$b->set_window_size($window);
	}
}

# Step3: for each chromosome, do statistics on windows and merge continuous windows into one region.
$time_tag = time;
my $global_data = {	# hash to hold global data.
	bedTr => \@tr_srbed,	# array of references to treatment SRBed objects.
	bedCo => \@co_srbed,	# array of references to control SRBed objects.
	bedBtr => \@btr_srbed,	# array of references to background treatment SRBed objects.
	bedBco => \@bco_srbed,	# array of references to background control SRBed objects.
	normTr => \@tr_norm,	# array of treatment normalization constants.
	normCo => \@co_norm,	# array of control normalization constants.
	normTr_ => \@tr_norm_,	# array of treatment normalization constants using #reads.
	normCo_ => \@co_norm_,	# array of control normalization constants using #reads.
	normBtr => \@btr_norm,	# array of background treatment normalization.
	normBco => \@bco_norm,	# array of background control normalization.
	bkgEnr => scalar @btr || scalar @bco,	# boolean tag for background enrichment.
	useBkg => $bkg,			# enrichment vs. background as filter.
	meanTr => \@tr_mean,	# mean bin count for treatment.
	meanCo => \@co_mean,	# mean bin count for control.
	stdTr => \@tr_std,		# bin count std for treatment.
	stdCo => \@co_std,		# bin count std for control.
	nsd => $nsd,			# number of std above mean as count cutoff.
	winSize => $window,		# window size.
	step => $step,			# moving step size.
	gap => $gap				# gap size.
};
my %chr_regList;	# hash table to hold region lists for each chromosome.
my %chr_winCov;	# hash table to hold window coverage for each chromosome.
if(!DEBUG) {
	# Create forkmanager for differential tests.
	$pm->run_on_finish(	# callback function to retrieve SRBed objects.
		sub{
			my($pid, $exit_code, $ident, $exit_signal, $core_dump, $dr) = @_;
			if(defined $dr){
				$chr_regList{$dr->{'chrom'}} = $dr->{'regList'};
				$chr_winCov{$dr->{'chrom'}} = $dr->{'winCov'};
			}else{
				# $chr_regList{$dr->{'chrom'}} = undef;
				# $chr_winCov{$dr->{'chrom'}} = undef;
				warn "Child process returned void data in differential tests ".
					"with exit code: $exit_code.\n";
			}
		}
	);
}
# Iterate through all chromosomes.
foreach my $r(@chrlen_ordered){
	my $slideWindow = new diffReps::SlideWindow;	# sliding window.
	my $regList = new diffReps::RegionList;	# list to store significant regions.
	my $winCov = new WindowCoverage($window, $step);	# class to calc window coverage.
	my $chrom = $r->{chrom};
	if(!DEBUG) { $pm->start and next; }
	# Create bin count and window count vectors for the current chromosome.
	# $window_num should be the same for all samples.
	my $window_num = create_chr_cntvec(\@tr_srbed, $chrom);
	if(!$window_num) {
		warn "Create count vector for treatment at: $chrom failed! Skip.\n";
		next ENDCHR;
	}
	if(!create_chr_cntvec(\@co_srbed, $chrom)) {
		warn "Create count vector for control at: $chrom failed! Skip.\n";
		next ENDCHR;
	}
	if(@btr_srbed) {
		if(!create_chr_cntvec(\@btr_srbed, $chrom)) {
			warn "Create count vector for background treatment at: $chrom failed! Skip.\n";
			next ENDCHR;
		}
	}
	if(@bco_srbed) {
		if(!create_chr_cntvec(\@bco_srbed, $chrom)) {
			warn "Create count vector for background control at: $chrom failed! Skip.\n";
			next ENDCHR;
		}
	}
	# Prepare window to slide.
	$slideWindow->set_chrom($global_data, $chrom, $window_num);
	my $signRegion = new diffReps::ChromaModSite;	# significant region.
	$winCov->init_movingCov();

	# Continue search for significant windows and expand the region.
	while($slideWindow->get_winN()>=0){
		if($slideWindow->isAboveCutoff($global_data)){
			$winCov->addcov();
			my $pval = $slideWindow->calc_stat($global_data, $meth, EPSILON);
			if($pval < $sign){	# Found a new significant window.
				if($signRegion->isWithinGap($global_data, $slideWindow)){	# if within gap, expand region.
					$signRegion->expand($global_data, $slideWindow);
				}
				else{	# No longer within gap: stop expanding; record old region and initialize a new one.
					if($signRegion->needStat()){
						#$signRegion->retrieve_norm_cnt($global_data);
						#$signRegion->retrieve_raw_sum($global_data);
						$signRegion->calc_stat($global_data, $meth, EPSILON);
					}
					# Nothing is done if current region is empty.
					$regList->add($signRegion);	# add region into list; clean current data structure.
					$signRegion->init($global_data, $slideWindow);	# initialize a region with current window.
				}
			}
		}
		$slideWindow->move_forward($global_data);
		$winCov->move();
	}
	# Print the last significant region left in memory.
	if($signRegion->needStat()){
		#$signRegion->retrieve_norm_cnt($global_data);
		#$signRegion->retrieve_raw_sum($global_data);
		$signRegion->calc_stat($global_data, $meth, EPSILON);
	}
	$regList->add_del($signRegion);
	
	ENDCHR: {
		# Delete bin count and window count vectors to save memory. 
		foreach my $r(@tr_srbed){
			$r->del_bin_count($chrom);
			$r->del_win_count($chrom);
		}
		foreach my $r(@co_srbed){
			$r->del_bin_count($chrom);
			$r->del_win_count($chrom);
		}
		if(@btr_srbed){
			foreach my $r(@btr_srbed){
				$r->del_bin_count($chrom);
				$r->del_win_count($chrom);
			}
		}
		if(@bco_srbed){
			foreach my $r(@bco_srbed){
				$r->del_bin_count($chrom);
				$r->del_win_count($chrom);
			}
		}
		if(DEBUG) {
			$chr_regList{$chrom} = $regList;
			$chr_winCov{$chrom} = $winCov;
		} else {
			my %rd = ('chrom' => $chrom, 'regList' => $regList, 'winCov' => $winCov);
			$pm->finish(0, \%rd);
		}
	}
}
if(!DEBUG) { $pm->wait_all_children; }
# Re-assemble data according to sorted chromosome names.
my $regList = new diffReps::RegionList;	# list to store significant regions.
my $winCov = new WindowCoverage($window, $step);	# class to calc window coverage.
foreach my $r(@chrlen_ordered){
	my $cn = $r->{'chrom'};
	if(exists $chr_regList{$cn} and exists $chr_winCov{$cn}) {
		$regList->append_list($chr_regList{$cn});
		$winCov->sumcov($chr_winCov{$cn});
	}
}
# Multiple testing correction.
my $N = floor($winCov->get_wincov() / $window);
$regList->adjPval($N);	# adjust P-value using estimated 'N'.
spit($hrep, '## FDR control information ##');
spit($hrep, '# Genome coverage for windows that passed cutoff', $winCov->get_wincov());
spit($hrep, '# N used in BH adjustment', $N);
spit($hrep, '# Detecting differential sites done in ' . (time - $time_tag) . ' seconds');
$hrep->flush;
# Dump all differential sites.
$regList->gen_header($hrep);	# result table header line.
$regList->output($global_data, $hrep);	# dump all found regions after P-value adjustment.
close $hrep;


# Step4: Annotate differential sites.
unless($noanno or $gname eq ''){
	`region_analysis.pl -i $report -r -d refseq -g $gname`;
}

# Step5: Look for hotspots.
unless($nohs){
	my $hotspot = $report . '.hotspot';
	`findHotspots.pl -d $report -o $hotspot`;
}

# Clean up - delete all temporary chromosome files.
if(@tr_srbed or @co_srbed or @btr_srbed or @bco_srbed) {
	foreach my $r( (@tr_srbed, @co_srbed, @btr_srbed, @bco_srbed) ) {
		$r->del_all_chrfiles();
	}
}

######### The END ##########


################ START Subroutines ###################

# Choose window and step sizes.
sub choose_winstep{
	my($mode,$win,$step) = @_;
	if($win <= 0 and $step <= 0){
		if($mode eq 'peak' or $mode eq 'p'){
			$win = 1000; $step = 100;
		}elsif($mode eq 'nucleosome' or $mode eq 'n'){
			$win = 200; $step = 20;
		}elsif($mode eq 'block' or $mode eq 'b'){
			$win = 10000; $step = 1000;
		}else{
			warn "Unspecified scanning mode. Use peak instead!\n";
			$win = 1000; $step = 100;
		}
	}elsif($step <= 0){
		if($win % 10 == 0){
			$step = $win / 10;
		}else{
			warn "Window size is not a multiple of 10. Set step size to be equal to window size.\n";
			$step = $win;
		}
	}else{
		if($win % $step != 0){
			warn "Window size is not a multiple of step size. Set step size to be equal to window size.\n";
			$step = $win;
		}
	}
	return ($win, $step);
}

# Get chrlen as an ordered array from chrlen hash table.
sub get_ordered_chrlen{
	my $ref_chrtbl = shift;
	my @a_chrlen;
	while(my($chrom, $chrlen) = each %{$ref_chrtbl}){
		push @a_chrlen, {'chrom' => $chrom, 'len' => $chrlen};
	}
	return sort {order_chr($a, $b)} @a_chrlen;
}

# Create bin and window count vectors for a chromosome.
# Return number of windows.
sub create_chr_cntvec{
	my($rbeds, $chrom) = @_;
	if(!@{$rbeds}) {return 0;}
	foreach my $b(@{$rbeds}){
		$b->chr_bin_count($chrom);
		$b->chr_win_count($chrom);
		return 0 unless defined $b->{chr_list}{$chrom}->{window_num};
	}
	return $rbeds->[0]->{chr_list}{$chrom}->{window_num};
}

# Go through a set of SRBed objects and calculate the mean, std and normalization constants.
sub mean_std_norm{
	my($robust,$rbeds,$chrlen) = @_;
	my $nsd = undef;	# cutoff num of std.
	if(@_ > 3) {$nsd = $_[3];}	# if defined, do norm after std cutoff.
	my $res = {	# result hash table.
		mean =>	[],	# array of means.
		std =>	[],	# array of stds.
		norm =>	[]	# array of normalization constants.
	};
	my @scales;	# array of scales.
	foreach my $b(@{$rbeds}) {
		push @scales, [];
	}
	#### Calculate mean and std. ####
	if($robust->{'std'}){	# standard estimation.
		my(@sums,@sumsqs);	# array of sum and sum of squares for all SRBed objects.
		foreach my $b(@{$rbeds}) {	# initialize arrays.
			push @sums, 0;
			push @sumsqs, 0;
		}
		my $nbin = 0;	# number of bins in the genome.
		foreach my $r(@{$chrlen}){
			my $chrom = $r->{chrom};
			foreach my $b(@{$rbeds}) {$b->chr_bin_count($chrom);}
			my $maxN = $rbeds->[0]->{chr_list}{$chrom}->{bin_num};	# num of bins in chrom.
			$nbin += $maxN;
			for my $i(0..$maxN-1){	# go through each bin.
				my @cnt;
				my $j = 0;	# iterator for SRBed objects.
				foreach my $b(@{$rbeds}) {
					push @cnt, $b->get_bin_count($chrom, $i);
					if($cnt[$j] > 0){
						$sums[$j] += $cnt[$j];
						$sumsqs[$j] += $cnt[$j]**2;
					}
					$j++;
				}
			}
		}
		foreach my $s(@sums) {push @{$res->{mean}}, $s/$nbin;}
		my $j = 0;	# iterator for SRBed objects.
		foreach my $s(@sumsqs) {push @{$res->{std}}, sqrt($s/$nbin - $res->{mean}[$j++]**2);}
	}else{	# robust estimation.
		my @sample_cnt;	# array of arrays of sampled counts.
		foreach my $b(@{$rbeds}) {	# initialize arrays.
			push @sample_cnt, [];
		}
		# Determine the random sampling rate.
		my $nbin = 0;	# number of bins in the genome.
		foreach my $r(@{$chrlen}){
			my $chrom = $r->{chrom};
			foreach my $b(@{$rbeds}) {$b->chr_bin_count($chrom);}
			my $maxN = $rbeds->[0]->{chr_list}{$chrom}->{bin_num};	# num of bins in chrom.
			$nbin += $maxN;
		}
		my $sam_rate = GENOM_SAMN / $nbin;
		# Take a sub-sample from whole genome ChIP-seq counts.
		foreach my $r(@{$chrlen}){
			my $chrom = $r->{chrom};
			my $maxN = $rbeds->[0]->{chr_list}{$chrom}->{bin_num};	# num of bins in chrom.
			for my $i(0..$maxN-1){	# go through each bin.
				my $rann = rand(1);
				my $j = 0;	# iterator for SRBed objects.
				foreach my $b(@{$rbeds}) {
					if($rann < $sam_rate){
						push @{$sample_cnt[$j++]}, $b->get_bin_count($chrom, $i);
					}
				}
			}
		}
		# Calcualte trimmed-mean and mad using sampled counts.
		foreach my $s(@sample_cnt){
			my $m = mean_r($s, 0, $robust->{'alpha'});
			push @{$res->{mean}}, $m;
			push @{$res->{std}}, mad($s, $m);
		}
	}

	#### Estimate norm contants. ####
	foreach my $r(@{$chrlen}){
		my $chrom = $r->{chrom};
		my $maxN = $rbeds->[0]->{chr_list}{$chrom}->{bin_num};	# num of bins in chrom.
		for my $i(0..$maxN-1){	# go through each bin.
			my @cnt;	# array to book-keep counts that pass cutoff.
			my $j = 0;	# iterator for SRBed objects.
			foreach my $b(@{$rbeds}) {
				my $c = $b->get_bin_count($chrom, $i);
				if(defined $nsd){
					last if $c < $res->{mean}[$j] + $nsd*$res->{std}[$j];
				}
				push @cnt, $c;
				$j++;
			}
			if(@cnt == @{$rbeds}){	# all counts pass std cutoff.
				my $gm = geomean(@cnt);	# bin geometric mean.
				if($gm > 0) { for $j(0..$#cnt) {push @{$scales[$j]}, $cnt[$j] / $gm;} }
			}
		}
		# Delete bin count vectors to free memory.
		foreach my $b(@{$rbeds}) {$b->del_bin_count($chrom);}
	}
	foreach my $s(@scales) {push @{$res->{norm}}, median(@{$s});}

	return $res;
}

# Split an array into two halves.
sub split_array{
	my($O,$len,$A,$B) = @_;
	@{$B} = @{$O};
	@{$A} = splice @{$B}, 0, $len;
}

# Dump variable values and descriptions to output.
sub spit{
	my $h = shift;	# output handle.
	my $desc = shift;
	if(@_ > 1){
		print $h "$desc: ", join("\t", @_), "\n";
	}
	elsif(@_==1){
		print $h "$desc: ", $_[0], "\n";
	}
	else{
		print $h "$desc\n";
	}
}

################## END Subroutines ####################


# A convenient class for calculating coverage of windows above count cutoff.
# Window coverage is finally used to calculate 'N' in BH P-value adjustment.
package WindowCoverage;
use strict;

sub new{
	die "Not enough arguments to set up window coverage class! Stop.\n" if @_ < 3;
	my($class,$winSize,$step) = @_;
	my $self = {
		movingCov => 0,
		windowCov => 0,
		winSize => $winSize,
		step => $step
	};
	bless $self, $class;
	return $self;
}

# Increment moving coverage counter when window moves.
sub move{
	my($self) = @_;
	# Add coverage if window moves. But cannot exceed window size.
	$self->{movingCov} += $self->{step} if $self->{movingCov} < $self->{winSize};
}

# Add window coverage for the first significant window on the chromosome.
sub init_movingCov{
	my $self = shift;
	$self->{movingCov} = $self->{winSize};
}

# Add accumulated coverage then reset counter.
sub addcov{
	my($self) = @_;
	$self->{windowCov} += $self->{movingCov};
	$self->{movingCov} = 0;
}

# Add the coverage from another obj into this one.
sub sumcov{
	my($self,$wc) = @_;
	$self->{'windowCov'} += $wc->{'windowCov'};
}

sub get_wincov{
	my $self = shift;
	return $self->{windowCov};
}








