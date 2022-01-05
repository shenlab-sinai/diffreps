Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg

NEWS
----
The diffReps project has migrated from Google code to Github. The original download files can be found in the **Downloads** folder.

NEWEST CHANGES & BUG FIXES
--------------------------

Please visit diffReps' discussion group for announcements and bug fixes:

<https://groups.google.com/forum/?fromgroups\#!forum/diffreps-discuss>

NEW REGION ANALYSIS PROGRAM
---------------------------

The region analysis program has spawn into an independent project that is now developed by Dr. Ningyi Shao. You may want to pay it a visit:

<https://github.com/shenlab-sinai/region_analysis>

It contains many more genomes than the default one in diffReps.

INTRODUCTION
------------

ChIP-seq is now widely used to profile the protein DNA interactions on a genome. It is of high interest to compare the differential enrichment of a histone mark or transcription factor between two contrasting conditions, such as disease vs. control. diffReps is developed to serve this purpose. It scans the whole genome using a sliding window, performing millions of statistical tests and report the significant hits. diffReps takes into account the biological variations within a group of samples and uses that information to enhance the statistical power. Considering biological variation is of high importance, especially for in vivo brain tissues.

In addition, diffReps has two easy-to-use tools. One is for quick genomic annotation of a differential site or peak list. Another is for finding chromatin modification hotspots (see below). diffReps is developed as an automated pipeline so that everything can be done in one command. In incorporates four different statistical tests, i.e. negative binomial, T-test, G-test and Chi-square test, for differential analysis. So no matter you have biological replicates or not, you can choose an appropriate test for your purpose.

I uploaded some slides about diffReps:

<div style="margin-bottom:5px">
<strong> <a href="http://www.slideshare.net/thefacultyl/diffreps-automated-chipseq-differential-analysis-package" title="diffReps: automated ChIP-seq differential analysis package" target="_blank">diffReps: automated ChIP-seq differential analysis package</a> </strong> from <strong><a href="http://www.slideshare.net/thefacultyl" target="_blank">Li Shen</a></strong>
</div>

PREREQUISITES
-------------

diffReps requires the following two CPAN modules:
```perl
Statistics::TTest
Math::CDF
Parallel::ForkManager
```
they can be downloaded and installed from CPAN. If you use cpanminus to install diffReps, they will be automatically installed.

Some systems have reported missing package `Time::HiRes`. If that is the case, it can be installed from CPAN or your choice of package manager. To test if it is already installed, use `perldoc Time::HiRes`.

INSTALLATION
------------

Installing diffReps is just like a standard PERL module. Basically you extract the package downloaded, go to the program directory and type the following commands:
```
perl Makefile.PL (Optional, PREFIX=your_perl_directory)
make
make test
make install
```
If you have root privileges, diffReps.pl will most likely be installed in `/usr/bin/`. If you specified PREFIX in Makefile, it will be installed in `your_perl_directory/bin`. Add `your_perl_directory/bin` to your PATH environmental variable, or copy diffReps.pl from `your_perl_directory/bin` to a directory that is already in PATH, such as `/home/yourname/bin`.

Alternative. If you have cpanminus installed, you can also install diffReps with one line command
```
cpanm diffReps-XXX.tar.gz
```
it will try to satisfy all the dependencies for you.

USE REGION ANALYSIS
-------------------

Please be notified that the included `refgene_getnearestgene` is built on Linux OS. If you are using Mac or Windows, you should download the corresponding program from cisgenome website: <http://www.biostat.jhsph.edu/~hji/cisgenome/> and over-write the one that is installed. To determine where the program is installed, use `which refgene_getnearestgene`.

The diffReps main program should run on all OS. If you do not want region analysis, you can turn it off in commands.

CHROMOSOME LENGTHS
------------------

It is important to supply diffReps with chromosome length information. diffReps requires that to bin the chromosomes into smaller sections. diffReps has a few genomes built-in so what you need to do is just give a genome name, such as mm9 or hg19. If the genome you are interested in is not already defined, you can give a text file for chorosome lengths. An example input is like
```
chr1    197195431
chr2    181748086
chr3    159599782
...
```

A NOTE ABOUT STATISTICAL TESTS
------------------------------

When you have biological replicates, Negative Binomial(NB) is the recommended test for differential analysis. An exact NB test is implemented in diffReps. Because NB distribution models discrete count data and over-dispersion among different samples, it appears to be an ideal model for ChIP-seq data. Many studies to date have used T-test on normalized counts for differential analysis. However, this is sub-optimal because normalized counts are NOT Normally distributed! As a result, detection power can be significantly degraded. Another caveat about T-test is that regions with very small counts may be picked up. Those regions should never pass cutoff because they don't have statistical significance. T-test ignores this fact because it simply treats them as continuous values. I still provide T-test in diffReps just for comparison purpose.

If your experiment doesn't contain biological replicates, you can choose between G-test and Chi-square test for differential analysis. They both give similar results but G-test is more recommended and has gained its popularity recently. See <http://en.wikipedia.org/wiki/G-test> for explanation. When they are chosen, diffReps performs a goodness-of-fit test on the normalized counts of treatment and control groups.

You can also use G-test or Chi-square test on data WITH biological replicates. An incentive of doing this is that this may give you more sensitivity but with a possibility of incurring false positives. diffReps automatically combines the biological replicates and generate a probablity vector accordingly. That means, if you have TWO replicates for treatment group and THREE replicates for control group, the probablity vector will be adjusted to reflect the replicate number difference.

### Using control/input samples

If you have control samples such as DNA inputs or IgG, you can give them to diffReps. It will calculate fold enrichment ratios for each differential site, which can be used for further filtering purposes. However, diffReps does NOT use those control samples for differential analysis.

DIFFERENTIAL SITES ANNOTATION
-----------------------------

diffReps includes a script for annotation of a differential sites list. By default, it will be evoked after diffReps finished running and annotate the differential sites based on their locations to the nearest genes. If no nearby genes can be found, it will also associate the differential sites with heterochromatic regions. A differential site will be assigned to one of the following categories:

| Region | Descriptions |
| ------ | ------------ |
| ProximalPromoter | +/- 250bp of TSS | 
| Promoter1k | +/- 1kbp of TSS | 
| Promoter3k | +/- 3kbp of TSS | 
| Genebody | Anywhere between a gene's promoter and up to 1kbp downstream of the TES. | 
| Genedeserts | Genomic regions that are depleted with genes and are at least 1Mbp long. | 
| Pericentromere | Between the boundary of a centromere and the closest gene minus 10kbp of that gene's regulatory region. | 
| Subtelomere | Similary defined as pericentromere. | 
| OtherIntergenic | Any region that does not belong to the above categories. |

The script can also be triggered manually. For example, if you want to annotate a differential list diff.h3k4me3.txt, you can use command like:
```
region_analysis.pl -i diff.h3k4me3.txt -r -d refseq -g mm9
```
will annotate the list using reference genome mm9 and the RefSeq database. The output will write to `diff.h3k4me3.txt.annotated`.

FINDING HISTONE MODIFICATION HOTSPOTS
-------------------------------------

The distance between two adjacent differential sites can be approximated by a Poisson distribution if they were positioned by random allocation. In reality, differential sites are often discovered to be spatially clustered together, forming so called chromatin modification hotspots. diffReps finds the hotspots by first building a null model on site-to-site distance, and then looking for regions that violate the null model with statistical significance using greedy search. diffReps reports the start, end positions as well as associated p-values and FDR of the hotspots. In addition, diffReps can accept more than one differential list of different histone marks as input, so that one can predict hotspots that show interaction between two or more histone marks.

By default, the finding hotspots routine will be called after diffReps finishes detecting differential sites. It will try to identify hotspots from the differential list just generated. The routine can also be used separately. For example, if you have two differential lists named diff.h3k4me3.txt and diff.polII.txt, you can look for hotspots that represent the interaction between H3K4me3 and Pol II using command like:
```
findHotspots.pl -d diff.h3k4me3.txt diff.polII.txt -o hotspot_k4.pol2.txt
```
will generate a hotspots list in `hotspot_k4.pol2.txt` file.

EXAMPLES
--------

diffReps requires input of BED files for ChIP-seq alignments for both treatment and control groups. BED files can be converted from any alignment format, such as BAM (Tip: you can use BedTools for this). An example of using diffReps for differential analysis is as follows
```
diffReps.pl -tr C1.bed C2.bed C3.bed -co S1.bed S2.bed S3.bed -gn mm9 \
-re diff.nb.txt -me nb
```
The output will be in diff.nb.txt file. If you want to specify your own chromosome/scaffold size in a text file such as chrlen.txt, you can use this: `-ch chrlen.txt OR --chrlen chrlen.txt`. By default, a sliding window of 1kbp is used with a moving step size of 100bp. There are other parameters that can be tuned for your data. Just type diffReps.pl in command console without specifying any arguments and hit Enter, you will see a usage summary.

RUNNING TIME
------------

The running time of diffReps totally depends on data and parameter settings. It can vary wildly between 30min and 10h. The most influential parameters on running time are window size and step size. The smaller the window size, the longer the running time which scales linearly.

HOW TO CITE
-----------

Shen L, Shao N-Y, Liu X, Maze I, Feng J, et al. (2013) diffReps: Detecting Differential Chromatin Modification Sites from ChIP-seq Data with Biological Replicates. PLoS ONE 8(6): e65598.

PAPERS THAT CITE DIFFREPS
-------------------------
A [list of papers](https://scholar.google.com/scholar?hl=en&as_sdt=5,33&cites=11143041154318133867&scipsc=&q=&scisbd=1) that cite diffReps found by Google Scholar.


CONTACT
-------

diffReps is developed by Dr. Li Shen at the Icahn School of Medicine at Mount Sinai. For technical issues, please post your questions to the discussion forum. If you want to collaborate with us, you may write to:  li.shen**AT**mssm.edu.

# TERMS OF USE
All data is free to use for non-commercial purposes. For commercial use please contact [MSIP](https://www.ip.mountsinai.org/).
