#!/usr/bin/perl -w
# Identify chromatin modification hotspots from differential analysis results.

use strict;
our $VERSION = '1.11.2';
use Getopt::Long;
use diffReps::DiffRes;

my $sign = 1e-3;
my $win1 = 1e6;
my $win2 = 5e6;
my $iformat = '11';
my @diffList;
my $output;
my $cmd = join(' ', $0, @ARGV);
if(@ARGV < 1 or !GetOptions('diff=s{1,}' => \@diffList,
							'output=s' => \$output,
							'pval:f' => \$sign,
							'win1:i' => \$win1,
							'win2:i' => \$win2,
							'infocol:s' => \$iformat)){
	print "Usage: findHotspots.pl --diff diff_file1..[diff_fileN] --output output_file \\\n";
	print "         [--pval p_value_cutoff] [--infocol info_cols_output]\n\n";
	print "Optional parameters(defaults in parentheses):\n";
	print "--pval(0.001)    P-value cutoff for significant clusters.\n";
	print "--win1(1000000)  1st window size for local lambda.\n";
	print "--win2(5000000)  2nd window size for local lambda.\n";
	print "--infocol(11)    Columns to extract as output. Use 0 to select none.\n";
	print "                 Default=11, which is the event column in diffReps output.\n";
	print "                 Format example: 4,6-8,11 will extract column 4,6,7,8 and 11.\n";
	exit 0;
}

if($win1 >= $win2) {
	print "1st window must be smaller than 2nd window.\n";
	exit 0;
}

## Start main procedure ##
my $diffRes = new diffReps::DiffRes($win1, $win2);
$diffRes->read_modsite(\@diffList, $iformat);	# load all differential lists.
$diffRes->sort_list();	# sort by genomic coordinates.
$diffRes->calc_glbL();	# global lambda(avg inter-distance between two sites).
if($diffRes->{'glbLambda'} == $diffRes->INF){
	die "Global lambda cannot be obtained. No hotspots need to be searched. Exit.\n";
}
my $pre_p = $diffRes->init_clust();	# initialize cluster by the 1st site.
my $clustList = new diffReps::ClustList();
$clustList->open_file($output);
$clustList->spit('# Command used: ' . $cmd);
$clustList->spit('# P-value cutoff: ' . $sign);
$clustList->print_header();
while(!$diffRes->is_list_end()){
	my $p = $diffRes->push_clust();
	if($p < $pre_p){
		$pre_p = $p;
		$p = $diffRes->shift_clust();	# remove the left most.
		while($p < $pre_p){	# iterative removing.
			$pre_p = $p;
			$p = $diffRes->shift_clust();
		}
		$pre_p = $diffRes->unshift_clust();	# last shift must be unsuccessful.
	}else{
		if($diffRes->pop_clust() < $sign){	# evaluate significance.
			$clustList->add($diffRes->get_clustinfo());
		}
		$pre_p = $diffRes->init_clust();
	}
}
if($diffRes->{'clust'}{'pval'} < $sign){	# don't forget the one left in list.
	$clustList->add($diffRes->get_clustinfo());
}
$clustList->adjPval($diffRes->get_siteN());
$clustList->output();
$clustList->done();	# clean up.




