#!/usr/bin/perl -w
#
# Analyze a set of regions by mapping them to various
# genomic features.
#

use strict;
use Getopt::Long;
require PJ::Database::Root;

# Program arguments.
sub help_msg{
	print "Usage: $0 --input in_reg_list [--rhead] [--help] [--database refseq(default)|ensembl] [--genome mm9(default)]\n";
	print "--input, -i      Input region file must assume the first 3 columns contain (chr,start,end)\n";
	print "\nOptional arguments:\n";
	print "--rhead, -r      Whether the input file contains column header\n";
	print "--help, -h       Print this help message\n";
	print "--database, -d   Choose database: refseq(default) or ensembl.\n";
	print "--genome, -g     Choose genome, currently supported: mm9(default), rn4, hg19.\n";
}

my $help = 0;
my $database = 'refseq';
my $genome = 'mm9';
my $head = 0;
my $regfile;
if(@ARGV < 1 or !GetOptions('help' => \$help,
							'input=s' => \$regfile,
							'database:s' => \$database,
							'genome:s' => \$genome,
							'rhead' => \$head)){
	&help_msg();
	exit 0;
}
if($help){
	&help_msg();
	exit 0;
}
die "Unrecognized database name.\n" unless $database eq 'refseq' or $database eq 'ensembl';
die "Unsupported genome.\n" unless $genome eq 'mm9' or $genome eq 'rn4' or $genome eq 'hg19';
my $dfile;
my $ggap_folder;
if(exists $INC{'PJ/Database/Root.pm'}){
	$ggap_folder = $INC{'PJ/Database/Root.pm'};
	$ggap_folder =~ s/Root.pm$//;
}else{
	die "Cannot locate Perl module PJ::Database::Root in search path. Genomic database must be installed to use region analysis.\n";
}
if($database eq 'refseq'){
	$dfile = $ggap_folder . $genome. '/' . $genome . '.RefSeq.refFlat.cisgenome.txt';
}elsif($database eq 'ensembl'){
	$dfile = $ggap_folder . $genome. '/' . $genome . '.ensembl.refFlat.cisgenome.txt';
}else{
}
die "Database file missing.\n" unless -f $dfile;

# Test whether the necessary programs are installed.
die "Cannot found program refgene_getnearestgene! Exit.\n" if system("which refgene_getnearestgene &> /dev/null");
my $ib_cmd;
if(!system("which intersectBed &> /dev/null")){
	$ib_cmd = "intersectBed";
}elsif(!system("which bedtools &> /dev/null")){
	$ib_cmd = "bedtools intersect";
}else{
	die "Cannot find intersectBed! BEDTools has not been installed.\n";
}

# Locate BED files for genedesert and pericentromere.
my $genedesert = $ggap_folder . $genome . '/' . $genome . '.genedeserts.bed';
my $pericentromere = $ggap_folder . $genome . '/' . $genome . '.pericentromeres.bed';
my $subtelomere = $ggap_folder . $genome . '/' . $genome . '.subtelomeres.bed';

# Remove comment lines from region file first.
my $regtmp = $regfile . '.tmp';
my $rm_comment = "grep -vP '\\W*#' $regfile > $regtmp";
`$rm_comment`;

# Convert regions to cod file for mapping.
if($head){
	`awk 'BEGIN{OFS="\t"} NR>1{print NR-1, \$1, \$2, \$3, "+"}' $regtmp > cod.tmp`;
	`awk 'BEGIN{OFS="\t"} NR>1{print \$1, \$2, \$3}' $regtmp > bed.tmp`;
} else{
	`awk 'BEGIN{OFS="\t"} {print NR, \$1, \$2, \$3, "+"}' $regtmp > cod.tmp`;
	`awk 'BEGIN{OFS="\t"} {print \$1, \$2, \$3}' $regtmp > bed.tmp`;
}

# Map regions to promoters and genebodies.
`refgene_getnearestgene -d $dfile -dt 1 -s mouse -i cod.tmp -o genes.tmp -r 0 -up 3000 -down 1000`;
# Map regions to genedeserts and pericentromeres.
`$ib_cmd -c -f 0.5 -a bed.tmp -b $genedesert > gd.tmp`;
`$ib_cmd -c -f 0.5 -a bed.tmp -b $pericentromere > pc.tmp`;
`$ib_cmd -c -f 0.5 -a bed.tmp -b $subtelomere > st.tmp`;
# Use the distance to TSS to designate promoter or genebody.
open GENES, '<', 'genes.tmp' or die "Itermediate output file missing: $!. Execution halt!\n";
open GD, '<', 'gd.tmp' or die "Itermediate output file missing: $!. Execution halt!\n";
open PC, '<', 'pc.tmp' or die "Itermediate output file missing: $!. Execution halt!\n";
open ST, '<', 'st.tmp' or die "Itermediate output file missing: $!. Execution halt!\n";
open ANNO, '>', 'anno.tmp' or die "Open temporary file for writing error: $!\n";
<GENES>;	# get rid of header.
print ANNO "GName\tTName\tStrand\tTSS\tTES\tFeature\tD2TSS\n" if $head;
while(<GENES>){
	chomp;
	# Read genedesert and pericentromere output.
	my $gdline = <GD>; chomp $gdline;
	my $pcline = <PC>; chomp $pcline;
	my $stline = <ST>; chomp $stline;
	# Parsing gene map results.
	my @cells = split /\t/;
	if(@cells < 11){	# not mapped to any gene.
		my @gdcells = split /\t/, $gdline;
		my @pccells = split /\t/, $pcline;
		my @stcells = split /\t/, $stline;
		if($gdcells[3] >= 1) {print ANNO "\t\t\t\t\tGenedesert\tNA\n";}
		elsif($pccells[3] >= 1) {print ANNO "\t\t\t\t\tPericentromere\tNA\n";}
		elsif($stcells[3] >= 1) {print ANNO "\t\t\t\t\tSubtelomere\tNA\n";}
		else {print ANNO "\t\t\t\t\tOtherIntergenic\tNA\n";}
		next;
	}
	my @output = ($cells[5], $cells[6], $cells[8], $cells[9], $cells[10]);
	# 3,4,9,10,11
	my $center2tss;
	if($cells[8] eq '+'){
		$center2tss = ($cells[2]+$cells[3])/2 - $cells[9];
	}else{
		$center2tss = $cells[10] - ($cells[2]+$cells[3])/2;
	}
	# NOTE: Proximal promoter is defined as a region +/-250bp of TSS.
	if(abs($center2tss) < 250){ 
		push @output, 'ProximalPromoter';
	} elsif(abs($center2tss) < 1000){
		push @output, 'Promoter1k';
	} elsif(abs($center2tss) < 3000){
		push @output, 'Promoter3k';
	} else{
		push @output, 'Genebody';
	}
	push @output, $center2tss;
	print ANNO join("\t", @output), "\n";
}
close GENES;
close ANNO;

my $annofile = $regfile . ".annotated";
`paste $regtmp anno.tmp > $annofile`;

# Clean up.
`rm -f anno.tmp bed.tmp cod.tmp gd.tmp genes.tmp pc.tmp st.tmp $regtmp`;
