#!/usr/bin/perl

###################################################################
# ExportVariants.pl
#
# Tool for exporting mutations from MLDBM to file
# 
#  Author: Giuseppe Narzisi 
#    Date: December 11, 2013
#
###################################################################

use warnings;
use strict;
#use POSIX;
use POSIX qw(strftime);
use FindBin qw($Bin);
use lib $Bin; # add $Bin directory to @INC
use Usage;
use Utils;
use SequenceIO;
use HashesIO;
use List::Util qw[min max];
use Getopt::Long;

use MLDBM::Sync;                       # this gets the default, SDBM_File
use MLDBM qw(DB_File Storable);        # use Storable for serializing
use MLDBM qw(MLDBM::Sync::SDBM_File);  # use extended SDBM_File, handles values > 1024 bytes
use Fcntl qw(:DEFAULT);                # import symbols O_CREAT & O_RDWR for use with DBMs

$|=1;

# required arguments
my $dbfile;
my $bedfile;

my $defaults = getDefaults();  

# defaults
my $_format  = $defaults->{format};
my $_SVtype  = $defaults->{SVtype};
my $_minCov  = $defaults->{min_cov};
my $_maxCov  = $defaults->{max_cov};
my $_maxchi2 = $defaults->{maxchi2};
my $_minCovRatio = $defaults->{covratio};
my $_intarget = 0;

# optional arguments (with default values)
my $format  = $_format;
my $SVtype  = $_SVtype;
my $minCov  = $_minCov;
my $maxCov  = $_maxCov;
my $intarget = $_intarget;
my $maxchi2 = $_maxchi2;
# min coverage ration (AleleCov/TotCov) for a mutation to be considered valid
# Mutation with coverage ration lower than this are considered sequencing errors.
my $minCovRatio = $_minCovRatio;

# Prob	chi2score
# 0.995	0.0000393
# 0.975	0.000982
# 0.20	1.642
# 0.10	2.706
# 0.05	3.841
# 0.025	5.024
# 0.02	5.412
# 0.01	6.635
# 0.005	7.879
# 0.002	9.550
# 0.001	10.828
#my $minchi2  = 10.828;

my %exons;
my %variants;

my $help;
my $VERBOSE = 0;

my $argcnt = scalar(@ARGV);
my $start_time = time;

#####################################################
#
# Command line options processing
#
GetOptions(
	'help!'    => \$help,
	'verbose!' => \$VERBOSE,
    
	# required parameters
    'db=s'  => \$dbfile,
    'bed=s' => \$bedfile,

	# optional paramters
    'format=s'   => \$format,
    'type=s'     => \$SVtype,
    'mincov=i'   => \$minCov,
    'maxcov=i'   => \$maxCov,
    'maxchi2=f'  => \$maxchi2,
    'covratio=f' => \$minCovRatio,
    'intarget!'  => \$intarget,

) or usageExport();

#####################################################
#
# Command line options processing
#
sub init()
{
	usageExport() if ($argcnt < 1);
	usageExport() if( $help );
	
	if( (!defined $dbfile) || (!defined $bedfile) ) { 
		print STDERR "Required parameter missing!\n";
		usageExport(); 
	}
	
	if ($minCov == 0) { $minCov = -1000000; }	
}

#####################################################

sub printParams {
	
	print STDERR "Parameters: \n";
	print STDERR "-- db-file: $dbfile\n";
	print STDERR "-- bed-file: $bedfile\n";
	print STDERR "-- output format: $format\n";
	print STDERR "-- SV type: $SVtype\n";
	print STDERR "-- minimum coverage: $minCov\n";
	print STDERR "-- maximum coverage: $maxCov\n";
	print STDERR "-- maximum chi-square score: $maxchi2\n";
	print STDERR "-- minimum coverage ratio: $minCovRatio\n";
	my $txt = "false";
	if($intarget) { $txt = "true"; }
	print STDERR "-- in target?: $txt\n\n";
}

## print de novo indels and check how how many are detected correctly
##########################################
sub printVariants {
	
	my $mode = $_[0];
		
	##  print header if file did not exist (created now)
	if($mode eq "annovar") { # annovar format
		print "#chr\tstart\tend\tref\tobs\tid\tsize\ttype\tavgKcov\tminKcov\tzygosity\taltKcov\tcovRatio\tchi2score\tinheritance\tbestState\tcovState\n"; 
	} 
	if($mode eq "vcf") { # vcf format
		print "##fileformat=VCFv4.1\n";
		my $date = strftime "%m/%d/%Y", localtime;
		print "##fileDate=$date\n";
		print "##source=scalpel$defaults->{version_num}\n";		
		#print "##reference=XXX\n";
				
		print "##INFO=<ID=AVGCOV,Number=1,Type=Float,Description=\"average k-mer coverage\">\n";
		print "##INFO=<ID=MINCOV,Number=1,Type=Integer,Description=\"minimum k-mer coverage of non-reference allele\">\n";
		print "##INFO=<ID=ALTCOV,Number=1,Type=Integer,Description=\"k-mer coverage of reference + any other allele (different from current non-reference) at locus\">\n";
		print "##INFO=<ID=ZYG,Number=1,Type=String,Description=\"zygosity\">\n";
		print "##INFO=<ID=COVRATIO,Number=1,Type=Float,Description=\"coverage ratio [(MINCOV)/(ALTCOV+MINCOV)]\">\n";
		print "##INFO=<ID=CHI2,Number=1,Type=Float,Description=\"chi-square score\">\n";
		print "##INFO=<ID=INH,Number=1,Type=String,Description=\"inheritance\">\n";
		print "##INFO=<ID=BESTSTATE,Number=1,Type=String,Description=\"state of the mutation\">\n";
		print "##INFO=<ID=COVSTATE,Number=1,Type=String,Description=\"coverage state of the mutation\">\n";
		print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
		print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"k-mer Depth\">\n";
		print "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"k-mer depth supporting reference/indel at the site\">\n";
		my $sample_name = "sample_name"; # replace by SM tag from @RG
		print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample_name\n";
	}
	elsif($mode eq "scalpel") { # scalpel format
		print "#ID\tchr\tpos\ttype\tlength\tavgKcov\tminKcov\tzygosity\tref\tobs\taltKcov\tloglikelihood\tchi2score\tinheritance\tbestState\tcovState\n"; 
	}
	
	my $num_snp = 0;
	my $num_ins = 0;
	my $num_del = 0;

	my $num_transitions = 0;
	my $num_transversions = 0;

	# prints denovo SVs to file
	foreach my $key (sort { bychrpos($a,$b,\%variants) } keys %variants) {
		#my $mut = $hash->{$key};
		my $mut = $variants{$key};
				
		my $family = $mut->{fam};
		$family = "na" if !defined $family;
		my $id = $mut->{id};
		$id = "na" if !defined $id;
		my $chr = $mut->{chr};
		my $pos = $mut->{pos};
		my $t = $mut->{type};
		my $l = $mut->{len};
		my $ref = $mut->{ref};
		my $qry = $mut->{seq}; 
		my $prevbpref = $mut->{prevbpref}; 
		my $prevbpalt = $mut->{prevbpalt}; 
		my $avgcov = $mut->{avgcov}; 
		my $mincov = $mut->{mincov}; 
		my $sta = $mut->{status};
		my $zyg = $mut->{zygosity};
		my $altcov = $mut->{altcov};
		my $inher = $mut->{inheritance};
		my $bestState = $mut->{bestState};
		$bestState = "na" if !defined $bestState;
		my $covState = $mut->{covState};
		$covState = "na" if !defined $covState;

		#**************** filtering ****************#		
		next if($mut->{mincov} < $minCov);
		next if($mut->{mincov} > $maxCov);
		if ($SVtype ne "all") {
			if($SVtype eq "indel") { next if($mut->{type} eq "snp"); }
			else { next if($mut->{type} ne $SVtype); }
		}
		#*******************************************#
				
		my $totcov = $altcov + $mincov;
		my $covRatio = sprintf('%.2f', ($mincov / $totcov) );
		
		#********************************************#
		next if($covRatio <= $minCovRatio); # skip sequencing errors
		#********************************************#

		# likelihood (-log p) of the mutation being error given the coverage at that locus		
		#my $loglikelihood = 0;
		#if( $prb != 1) { $loglikelihood = sprintf('%.2f', -10*log10(1.0 - $covRatio)); }
		
		#Chi-squared Test
		my $o1 = $mincov; # observed 1
		my $o2 = $altcov; # observed 2
		my $e1 = $totcov/2; # expected 1
		my $e2 = $totcov/2; # expected 2
		
		my $term1 = (($o1-$e1)*($o1-$e1))/$e1;
		my $term2 = (($o2-$e2)*($o2-$e2))/$e2;
		
		my $chi2Score = sprintf('%.2f', $term1 + $term2);		
		if($zyg eq "hom") { $chi2Score = 0; }
		
		#********************************************#
		next if($chi2Score > $maxchi2);
		#********************************************#
		
		my $annovar_ref = $ref;
		my $annovar_qry = $qry;
		my $vcf_ref = $prevbpref . $ref;
		my $vcf_qry = $prevbpalt . $qry;
		
		if($sta eq "ok") { ## only report clean indels...
			
			if($t eq "snp") { 
				$num_snp++;
				
				#transitions: A-G, C-T
				if ( ($ref eq "A") && ($qry eq "G") ||
				 	 ($ref eq "G") && ($qry eq "A") ||
					 ($ref eq "C") && ($qry eq "T") ||
					 ($ref eq "T") && ($qry eq "C") ) {
						$num_transitions++;
					}
				else { $num_transversions++; }
					
			}
			if($t eq "ins") { 
				$num_ins++; 
				$annovar_ref = "-";
				$vcf_ref = $prevbpref;
			}
			if($t eq "del") { 
				$num_del++; 
				$annovar_qry = "-";
				$vcf_qry = $prevbpalt;
			}
			
			my $start = $pos;
			my $end = $pos;
			if($t eq "del") { $end = $start+$l-1; }
			my $str;
			if($mode eq "annovar") { # annovar format
				$str = sprintf("%s\t%d\t%d\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $chr, $start, $end, $annovar_ref, $annovar_qry, $id, $l, $t, $avgcov, $mincov, $zyg, $altcov, $covRatio, $chi2Score, $inher, $bestState, $covState);
			}
			elsif($mode eq "vcf") { # vcf format	
				#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
 				my $format_str = "GT:AD:DP";
  				my $gt = $zyg eq "het" ? "0/1" : "1/1";
  				my $format_val = "$gt:$altcov,$mincov:$totcov";
				$covState=~s/ /,/g; # no spaces allowed by vcf spec
  				my $info = "AVGCOV=$avgcov;MINCOV=$mincov;ALTCOV=$altcov;ZYG=$zyg;COVRATIO=$covRatio;CHI2=$chi2Score;INH=$inher;BESTSTATE=$bestState;COVSTATE=$covState";
  				if($t eq "snp") { $str = sprintf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n", $chr, $start, ".", $ref, $qry, ".", "PASS", $info); }
  				else            { $str = sprintf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $chr, $start-1, ".", $vcf_ref, $vcf_qry, ".", "PASS", $info, $format_str, $format_val); }
			}
			else { 
				$str = sprintf("%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $chr, $pos, $t, $l, $avgcov, $mincov, $zyg, $ref, $qry, $id, $altcov, $covRatio, $chi2Score, $inher, $bestState, $covState);
			}
			print  "$str";
		}
 	}	
	
	my $ti_tv_ratio = 0;
	if ($num_transversions != 0) { $ti_tv_ratio = $num_transitions/$num_transversions; }
	print STDERR "Ti/Tv ratio: $ti_tv_ratio\n" if($VERBOSE);

	my $num_valid = $num_ins+$num_del+$num_snp;
	print STDERR "Total number of mutations: $num_valid\n" if($VERBOSE);
	
	#if($VERBOSE) {
		print STDERR "[#SNPs: $num_snp | #Ins: $num_ins | #Del: $num_del | Tot: $num_valid]\n";
	#}
}

## do the job
##########################################

init();
printParams() if ($VERBOSE);
loadExonsBed("$bedfile", \%exons, 0, $VERBOSE);
loadDB("$dbfile", \%variants, \%exons, $intarget);
printVariants("$format");

##########################################

my $time_taken = time - $start_time;

#if($VERBOSE) {
	elapsedTime($time_taken, "ExportVariants");
#}
