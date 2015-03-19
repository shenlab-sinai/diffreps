use 5.006;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
my @scripts = qw( diffReps.pl elandsorted2bed.pl export2sorted.pl rmdup_sorted_bed.pl sam2bed.pl findHotspots.pl refgene_getnearestgene region_analysis.pl );
WriteMakefile(
    NAME              => 'diffReps',
    VERSION_FROM      => 'bin/diffReps.pl', # finds $VERSION
  	EXE_FILES         => [ map {'bin/' . $_} @scripts ],
    PREREQ_PM         => {
  		'Statistics::TTest' => 0,
  		'Math::CDF' => 0,
      'Parallel::ForkManager' => 0,
      'Time::HiRes' => 0
  	}, # e.g., Module::Name => 1.1
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (#ABSTRACT_FROM  => 'lib/diffReps.pm', # retrieve abstract from module
  	   ABSTRACT       => 'Detecting differential chromatin modification sites from ChIP-seq data',
       AUTHOR         => 'Li Shen <li.shen@mssm.edu>') : ()
    ),
);
