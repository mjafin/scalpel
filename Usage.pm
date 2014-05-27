package Usage;

###################################################################
# Usage
#
# Package with usage description of each operation mode 
# 
#  Author: Giuseppe Narzisi 
#    Date: December 11, 2013
#
###################################################################

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(getDefaults usageDenovo usageSingle usageSomatic usageExport);

use strict;
use warnings;
use POSIX;

my %defaults = (
	kmer 				=> 25,
	cov_threshold 		=> 5,
	tip_cov_threshold	=> 2,
	covratio			=> 0.01,
	radius	 			=> 100,
	windowSize 			=> 400,
	delta				=> 100,
	map_qual			=> 1,
	maxmismatch 		=> 3,
	min_cov 			=> 5,
	max_cov				=> 1000000,
	maxchi2				=> 1000000,
	outratio			=> 0.05,
	WORK 				=> "./outdir",
	MAX_PROCESSES 		=> 1,
	sample		 		=> "ALL",
	selected 			=> "null",
	pathlimit 			=> 100000,
	format				=> "vcf",
	SVtype				=> "all",
	version_num			=> "0.1.1(beta)",
);

#####################################################
#
# export defaults
#
sub getDefaults { return \%defaults; }


#####################################################
#
# Message about the program and how to use it
#
sub usageSingle {	

my $name = $_[0];
print STDERR <<END;

usage: $name --bam <BAM file> --bed <BED file> --ref <FASTA file> [OPTIONS]
	
Detect INDELs in one single dataset (e.g., one individual).
	
OPTIONS:

    --help             : this (help) message
    --verbose          : verbose mode

  Required:
    --bam <BAM file>   : BAM file with the reference-aligned reads
    --bed <BED file>   : BED file with list of exome-target coordinates
    --ref <FASTA file> : reference genome in FASTA format

  Optional:
    --kmer <int>       : k-mer size [default $defaults{kmer}]
    --covthr <int>     : threshold used to select source and sink [default $defaults{cov_threshold}]
    --lowcov <int>     : threshold used to remove low-coverage nodes [default $defaults{tip_cov_threshold}]
    --covratio <float> : minimum coverage ratio for sequencing errors (default: $defaults{covratio})
    --radius <int>     : left and right extension (in base-pairs) [default $defaults{radius}]
    --window <int>     : window-size of the region to assemble (in base-pairs) [default $defaults{windowSize}]
    --step <int>       : delta shift for the sliding window (in base-pairs) [default$defaults{delta}]
    --mapscore <int>   : minimum mapping quality for selecting reads to assemble [default $defaults{map_qual}]
    --pathlimit <int>  : limit number of sequence paths [default $defaults{pathlimit}]
    --mismatches <int> : max number of mismatches in near-perfect repeat detection [default $defaults{maxmismatch}]
    --dir <directory>  : output directory [default $defaults{WORK}]
    --numprocs <int>   : number of parallel jobs (1 for no parallelization) [default $defaults{MAX_PROCESSES}]
    --sample <string>  : only process reads/fragments in sample [default $defaults{sample}]
    --coords <file>    : file with list of selected locations to examine [default $defaults{selected}]

  Output:
    --format           : export mutations in selected format (annovar | vcf) [default $defaults{format}]
    --intarget         : export mutations only inside the target regions from the BED file
    --mincov <int>     : minimum coverage for exporting mutation to file [default $defaults{min_cov}]
    --outratio <float> : minimum coverage ratio for exporting mutation to file (default: $defaults{outratio})

  Note 1: the list of detected INDELs is saved in file: OUTDIR/variants.*.indel.*
  where OUTDIR is the output directory selected with option "--dir" [default $defaults{WORK}]

  Note 2: the input reference file (option "--ref") must be the same one that was used to create the BAM file.

END
exit;
}

#####################################################
#
# Message about the program and how to use it
#
sub usageDenovo {
	
my $name = $_[0];
print STDERR <<END;

usage: $name --dad <BAM file> --mom <BAM file> --aff <BAM file> --sib <BAM file> --bed <BED file> --ref <FASTA file> [OPTIONS]

Detect de novo INDELs in a family of four individuals (mom, dad, aff, sib).

OPTIONS:

    --help             : this (help) message
    --verbose          : verbose mode

  Required:
    --dad <BAM file>   : father BAM file
    --mom <BAM file>   : mother BAM file
    --aff <BAM file>   : affected child BAM file
    --sib <BAM file>   : sibling BAM file
    --bed <BED file>   : BED file with list of exome-target coordinates
    --ref <FASTA file> : reference genome in FASTA format

  Optional:
    --kmer <int>       : k-mer size [default $defaults{kmer}]
    --covthr <int>     : threshold used to select source and sink [default $defaults{cov_threshold}]
    --lowcov <int>     : threshold used to remove low-coverage nodes [default $defaults{tip_cov_threshold}]
    --covratio <float> : minimum coverage ratio for sequencing errors (default: $defaults{covratio})
    --radius <int>     : left and right extension (in base-pairs) [default $defaults{radius}]
    --window <int>     : window-size of the region to assemble (in base-pairs) [default $defaults{windowSize}]
    --step <int>       : delta shift for the sliding window (in base-pairs) [default $defaults{delta}]
    --mapscore <int>   : minimum mapping quality for selecting reads to assemble [default $defaults{map_qual}]
    --mismatches <int> : max number of mismatches in near-perfect repeat detection [default $defaults{maxmismatch}]
    --dir <directory>  : output directory [default $defaults{WORK}]
    --numprocs <int>   : number of parallel jobs (1 for no parallelization) [default $defaults{MAX_PROCESSES}]
    --coords <file>    : file with list of selected coordinates to examine [default $defaults{selected}]

  Output:
    --format           : export mutations in selected format (annovar | vcf) [default $defaults{format}]
    --intarget         : export mutations only inside the target regions from the BED file
    --mincov <int>     : minimum coverage for exporting mutation to file [default $defaults{min_cov}]
    --outratio <float> : minimum coverage ratio for exporting mutation to file (default: $defaults{outratio})

  Note 1: the list of de novo INDELs is saved in file: OUTDIR/denovos.*.indel.*
  where OUTDIR is the output directory selected with option "--dir" [default $defaults{WORK}]

  Note 2: the input reference file (option "--ref") must be the same one that was used to create the BAM file.

END
exit;
}

#####################################################
#
# Message about the program and how to use it
#
sub usageSomatic {

my $name = $_[0];
print STDERR <<END;

usage: $name --normal <BAM file> --tumor <BAM file> --bed <BED file> --ref <FASTA file> [OPTIONS]

Detect somatic INDELs in a tumor/sample pair

OPTIONS:

    --help                : this (help) message
    --verbose             : verbose mode

  Required:
    --normal <BAM file>   : normal BAM file
    --tumor  <BAM file>   : tumor BAM file
    --bed    <BED file>   : BED file with list of exome-target coordinates
    --ref    <FASTA file> : reference genome in FASTA format

  Optional:
    --kmer <int>          : k-mer size [default $defaults{kmer}]
    --covthr <int>        : threshold used to select source and sink [default $defaults{cov_threshold}]
    --lowcov <int>        : threshold used to remove low-coverage nodes [default $defaults{tip_cov_threshold}]
    --covratio <float>    : minimum coverage ratio for sequencing errors (default: $defaults{covratio})
    --radius <int>        : left and right extension (in base-pairs) [default $defaults{radius}]
    --window <int>        : window-size of the region to assemble (in base-pairs) [default $defaults{windowSize}]
    --step <int>          : delta shift for the sliding window (in base-pairs) [default $defaults{delta}]
    --mapscore <int>      : minimum mapping quality for selecting reads to assemble [default $defaults{map_qual}]
    --mismatches <int>    : max number of mismatches in near-perfect repeat detection [default $defaults{maxmismatch}]
    --dir <directory>     : output directory [default $defaults{WORK}]
    --numprocs <int>      : number of parallel jobs (1 for no parallelization) [default $defaults{MAX_PROCESSES}]
    --coords <file>       : file with list of selected coordinates to examine [default $defaults{selected}]

  Output:
    --format              : export mutations in selected format (annovar | vcf) [default $defaults{format}]
    --intarget            : export mutations only inside the target regions from the BED file
    --mincov <int>        : minimum coverage for exporting mutation to file [default $defaults{min_cov}]
    --outratio <float>    : minimum coverage ratio for exporting mutation to file (default: $defaults{outratio})

  Note 1: the list of somatic INDELs is saved in file: OUTDIR/somatic.*.indel.* 
  where OUTDIR is the output directory selected with option "--dir" [default $defaults{WORK}]

  Note 2: the input reference file (option "--ref") must be the same one that was used to create the BAM file.

END
exit;
}


#####################################################
#
# Message about this program and how to use it
#
sub usageExport() {
print STDERR <<END;

usage: ExportVariants.pl --db <file> --bed <BED file> [OPTIONS]

OPTIONS:

    --help             : this (help) message
    --verbose          : verbose mode

  Required:
    --db <file>        : Database of mutations
    --bed <file>       : BED file with list of exome-target coordinates
  
  Optional:
    --format <text>    : output format for variants (annovar | vcf) [default $defaults{format}]
    --type <text>      : mutation type (snp, del, ins, indel, all: everything) [default $defaults{SVtype}]
    --mincov <int>     : minimum coverage for a mutation to be exported  [default $defaults{min_cov}]
    --maxcov <int>     : maximum coverage for a mutation to be exported  [default $defaults{max_cov}]
    --maxchi2 <float>  : maximum chi-square score for a mutation to be exported [default $defaults{maxchi2}]
    --covratio <float> : minimum coverage ratio (AlleleCov/TotCov) for a mutation to be exported to file [default $defaults{outratio}]
    --intarget         : export mutations only inside the target regions from the BED file

  Supported output formats:
    1. annovar
    2. vcf

END
exit;
}

1;