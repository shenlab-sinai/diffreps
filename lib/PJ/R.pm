package PJ::R;
our $VERSION = '0.01';

use strict;

our $R_calc_cov = '
	config.covCalc <- list(folder=".", pattern="INPUT", output="OUTPUT", adj=1.0)
	library(ShortRead)
	library(BSgenome)
	read.raw <- sapply(list.files(config.covCalc$folder, config.covCalc$pattern),
	  function(filename) readAligned(config.covCalc$folder, filename, type="TYPE"),
	  simplify=FALSE)
	genome <- "GENOME"
	load("CHRLEN")
	fragLen <- FRAG
	sumUpCoverage <- function( lanes, seqLens, fragmentLength )
	{
	  res <- NULL
	  for( i in 1:length(lanes) ) {
		filteredReads <- lanes[[i]][chromosome(lanes[[i]]) %in% names(seqLens)]
		# Calculate coverage for this lane.
		cvg <- coverage( filteredReads, width = seqLens,
		  extend = as.integer(fragmentLength) - width(filteredReads) )
		if( is.null( res ) )	# is the 1st lane?
		  res <- cvg
		else {	# if not, add this lane to existing results.
		  stopifnot( all( names(res) == names(cvg) ) )
		  for( seq in names(res) )
			res[[seq]] <- res[[seq]] + cvg[[seq]]
		}
	  }
	  res
	}
	read.coverage = sumUpCoverage(read.raw, seqlens, fragLen)
	nreads <- sum(sapply(read.raw, length))
	read.coverage.n <- GenomeData(lapply(read.coverage, function(r) r/nreads*1000000*config.covCalc$adj))
	save(genome, nreads, read.coverage, read.coverage.n, file=config.covCalc$output)';



1;
