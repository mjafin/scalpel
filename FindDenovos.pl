#!/usr/bin/perl

###################################################################
# FindDenovos.pl
#
# Tool for detecting de novo mutations in one family 
# using microassembly
# 
#  Author: Giuseppe Narzisi 
#    Date: December 11, 2013
#
###################################################################

use warnings;
use strict;
use POSIX;
use FindBin qw($Bin);
use lib $Bin; # add $Bin directory to @INC
use Usage;
use Utils qw(:DEFAULT $findVariants $findDenovos $findSomatic $exportTool $bamtools);
use HashesIO;
use Parallel::ForkManager;
use List::Util qw[min max];
#use Math::Random qw(:all);
#use Sys::CPU;
use Getopt::Long;
use File::Spec;

#use Tie::Cache;
use MLDBM::Sync;                       # this gets the default, SDBM_File
use MLDBM qw(DB_File Storable);        # use Storable for serializing
use MLDBM qw(MLDBM::Sync::SDBM_File);  # use extended SDBM_File, handles values > 1024 bytes
use Fcntl qw(:DEFAULT);                # import symbols O_CREAT & O_RDWR for use with DBMs

$|=1;

# required arguments
my $BAMDAD;
my $BAMMOM;
my $BAMAFF;
my $BAMSIB;
my $BEDFILE;
my $REF;

my $defaults = getDefaults();  

# optional arguments
my $kmer = $defaults->{kmer};
my $cov_threshold = $defaults->{cov_threshold};
my $tip_cov_threshold = $defaults->{tip_cov_threshold};
my $radius = $defaults->{radius};
my $windowSize = $defaults->{windowSize};
my $delta = $defaults->{delta};
my $map_qual = $defaults->{map_qual};
my $maxmismatch = $defaults->{maxmismatch};
my $min_cov = $defaults->{min_cov};
my $covratio = $defaults->{covratio};
my $outratio = $defaults->{outratio};
my $WORK  = $defaults->{WORK};
my $MAX_PROCESSES = $defaults->{MAX_PROCESSES};
my $selected = $defaults->{selected};
my $outformat = $defaults->{format};
my $intarget;

my $help = 0;
my $VERBOSE = 0;

my %candidates;
my %exons;
my %locations;
my %families;

my %fatherCov;
my %motherCov;
my %sibCov;
my %selfCov;

my %denovoL2K;
my %inheritedL2K;
my %fatherL2K;
my %motherL2K;
my %sibL2K;
my %selfL2K;

my %fatherSVs;
my %motherSVs;
my %sibSVs;
my %selfSVs;
my %denovoSVs;
my %inheritedSVs;

my @members = qw/mom dad aff sib/;
my @kids = qw/aff sib/;

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
    'dad=s' => \$BAMDAD,
    'mom=s' => \$BAMMOM,
    'aff=s' => \$BAMAFF,
    'sib=s' => \$BAMSIB,
    'bed=s' => \$BEDFILE,
    'ref=s' => \$REF,

	# optional paramters
    'kmer=i'       => \$kmer,
    'covthr=i'     => \$cov_threshold,
    'lowcov=i'     => \$tip_cov_threshold,
    'mincov=i'     => \$min_cov,
    'covratio=f'   => \$covratio,
    'radius=i'     => \$radius,
    'window=i'     => \$windowSize,
    'step=i'       => \$delta,
    'mapscore=i'   => \$map_qual,
    'mismatches=i' => \$maxmismatch,
    'dir=s'        => \$WORK,
    'numprocs=i'   => \$MAX_PROCESSES,
    'coords=s'     => \$selected,
    'format=s'     => \$outformat,

	# ouptut parameters
	'intarget!'    => \$intarget,
	'outratio=f'   => \$outratio,

) or usageDenovo("FindDenovo.pl");

#####################################################
#
# Initilization 
#
sub init()
{	
	usageDenovo("FindDenovo.pl") if ($argcnt < 1);
	usageDenovo("FindDenovo.pl") if( $help );
	
	if( (!defined $BAMDAD) || (!defined $BAMMOM) || (!defined $BAMAFF) || (!defined $BAMSIB) || (!defined $BEDFILE) || (!defined $REF) ) { 
		print STDERR "Required parameter missing!\n";
		usageDenovo("FindDenovo.pl"); 
	}
	
	mkdir $WORK if ! -r $WORK;
	
	if($MAX_PROCESSES == 1) { $MAX_PROCESSES = 0; }
}

#####################################################

sub printParams { 

	my $file = $_[0];
	
	print STDERR "-- Print paramters to $file\n";
	
	open PFILE, "> $file" or die "Can't open $file ($!)\n";
	
	print PFILE "Parameters: \n";
	print PFILE "-- K-mer size (bp): $kmer\n";
	print PFILE "-- coverage threshold: $cov_threshold\n";
	print PFILE "-- low coverage threshold: $tip_cov_threshold\n";
	print PFILE "-- size (bp) of the left and right extensions (radius): $radius\n";
	print PFILE "-- window-size size (bp): $windowSize\n";
	print PFILE "-- step size (bp): $delta\n";
	print PFILE "-- minimum mapping-quality: $map_qual\n";
	print PFILE "-- minimum coverage for denovo mutation: $min_cov\n";
	print PFILE "-- minimum coverage ratio for sequencing error: $covratio\n";
	print PFILE "-- minimum coverage ratio for exporting mutation: $outratio\n";
	print PFILE "-- max number of mismatches for near-perfect repeats: $maxmismatch\n";
	print PFILE "-- output directory: $WORK\n";
	print PFILE "-- max number of parallel jobs: $MAX_PROCESSES\n";
	print PFILE "-- father BAM file: $BAMDAD\n";
	print PFILE "-- mother BAM file: $BAMMOM\n";
	print PFILE "-- affected BAM file: $BAMAFF\n";
	print PFILE "-- sibling BAM file: $BAMSIB\n";
	print PFILE "-- bedfile: $BEDFILE\n";
	print PFILE "-- reference file: $REF\n";
	print PFILE "-- file of selected coordinates: $selected\n";
	print PFILE "-- output format for variants: $outformat\n";
	if($intarget) { print PFILE "-- output variants in target? yes\n"; }
	else { print PFILE "-- output variants in target? no\n"; }
	
	close PFILE;
}

## load mutations and coverage info
#####################################################

sub loadHashes {
	
	print STDERR "-- Loading mutations and coverage info\n";
		
	# clear hashes 
	for my $k (keys %fatherSVs) { delete $fatherSVs{$k}; }
	for my $k (keys %motherSVs) { delete $motherSVs{$k}; }
	for my $k (keys %sibSVs)    { delete $sibSVs{$k}; }
	for my $k (keys %selfSVs)   { delete $selfSVs{$k}; }
	
	for my $k (keys %fatherCov) { delete $fatherCov{$k}; }
	for my $k (keys %motherCov) { delete $motherCov{$k}; }
	for my $k (keys %selfCov) { delete $selfCov{$k}; }
	for my $k (keys %sibCov) { delete $sibCov{$k}; }
	
	for my $k (keys %fatherL2K) { delete $fatherL2K{$k}; }
	for my $k (keys %motherL2K) { delete $motherL2K{$k}; }
	for my $k (keys %sibL2K)    { delete $sibL2K{$k}; }
	for my $k (keys %selfL2K)   { delete $selfL2K{$k}; }
	
	for my $k (keys %denovoL2K) { delete $denovoL2K{$k}; }
	for my $k (keys %inheritedL2K) { delete $inheritedL2K{$k}; }
	
	# load databases of mutations
	foreach my $ID (@members) {
		print STDERR "load variants for $ID...\n";
		if($ID eq "dad") { loadDB("$WORK/$ID/variants.db", \%fatherSVs, \%exons, 0); }
		if($ID eq "mom") { loadDB("$WORK/$ID/variants.db", \%motherSVs, \%exons, 0); }
		if($ID eq "sib") { loadDB("$WORK/$ID/variants.db", \%sibSVs, \%exons, 0); }
		if($ID eq "aff") { loadDB("$WORK/$ID/variants.db", \%selfSVs, \%exons, 0); }
	}
	
	# load reference coverage info
	foreach my $ID (@members) {	
		print STDERR "load coverage for $ID...\n";
		if($ID eq "dad") { loadCov("$WORK/$ID/refcov.txt", \%fatherCov); }
		if($ID eq "mom") { loadCov("$WORK/$ID/refcov.txt", \%motherCov); }
		if($ID eq "sib") { loadCov("$WORK/$ID/refcov.txt", \%sibCov); }
		if($ID eq "aff") { loadCov("$WORK/$ID/refcov.txt", \%selfCov); }
	}
	
	# load loc2key info
	foreach my $ID (@members) {	
		print STDERR "load loc2key for $ID...\n";
		if($ID eq "dad") { loadLoc2Key("$WORK/$ID/loc2keys.txt.gz", \%fatherL2K); }
		if($ID eq "mom") { loadLoc2Key("$WORK/$ID/loc2keys.txt.gz", \%motherL2K); }
		if($ID eq "sib") { loadLoc2Key("$WORK/$ID/loc2keys.txt.gz", \%sibL2K); }
		if($ID eq "aff") { loadLoc2Key("$WORK/$ID/loc2keys.txt.gz", \%selfL2K); }
	}
	
	#remove reference coverage files:
	foreach my $ID (@members) {	
		print STDERR "remove coverage file for $ID...\n";
		runCmd("remove coverage file", "rm $WORK/$ID/refcov.txt.gz");
	}
	
}

## call variants on each family memeber
#####################################################
sub callSVs {
	
	print STDERR "-- Detect mutations on each family member\n";
		
	foreach my $ID (@members) {

		print STDERR "-- Processing $ID\n";
		
		my $bam_abs_path;
		if($ID eq "dad") { $bam_abs_path = File::Spec->rel2abs($BAMDAD); }
		if($ID eq "mom") { $bam_abs_path = File::Spec->rel2abs($BAMMOM); }
		if($ID eq "aff") { $bam_abs_path = File::Spec->rel2abs($BAMAFF); }
		if($ID eq "sib") { $bam_abs_path = File::Spec->rel2abs($BAMSIB); }
		
		my $command = "$findVariants ".
			"--bam $bam_abs_path ".
			"--bed $BEDFILE ".
			"--ref $REF ".
			"--kmer $kmer ". 
			"--covthr $cov_threshold ". 
			"--lowcov $tip_cov_threshold ".
			"--covratio $covratio ".
			"--radius $radius ". 
			"--window $windowSize ". 
			"--step $delta ".
			"--mapscore $map_qual ".
			"--mismatches $maxmismatch ".
			"--dir $WORK/$ID ".
			"--sample ALL ".
			"--numprocs $MAX_PROCESSES ".
			"--coords $selected ".
			"--format $outformat ".
			"--cov2file";
		if($intarget) { $command .= " --intarget"; }	
		if($VERBOSE) { $command .= " --verbose"; }
	
		print STDERR "Command: $command\n" if($VERBOSE);
		runCmd("findVariants", "$command");
	}
}

## find denovo events
#####################################################
sub findDenovos {
	
	print STDERR "-- Finding denovo mutations\n";
		
	if (-e "$WORK/denovos.db.dir") { runCmd("remove database files", "rm $WORK/denovos.db.*"); }
	if (-e "$WORK/inherited.db.dir") { runCmd("remove database files", "rm $WORK/inherited.db.*"); }
	
	my $denovo_dbm_obj = tie %denovoSVs, 'MLDBM::Sync', "$WORK/denovos.db", O_CREAT|O_RDWR, 0640; 
	my $inherited_dbm_obj = tie %inheritedSVs, 'MLDBM::Sync', "$WORK/inherited.db", O_CREAT|O_RDWR, 0640; 
	
	$denovo_dbm_obj->Lock;
	$inherited_dbm_obj->Lock;
	
	my $vhash;
	my $l2k;
	foreach my $ID (@kids) {
		
		if($ID eq "sib") { $vhash = \%sibSVs; $l2k = \%sibL2K; }
		if($ID eq "aff") { $vhash = \%selfSVs; $l2k = \%selfL2K; }
		
		foreach my $k (keys %$vhash) {
			my $mut = $vhash->{$k};
			
			# parse bestState and update info
			parseBestStateQuad($mut);
			
			my $chr = $mut->{chr};
			my $pos = $mut->{pos};
			my $key = "$chr:$pos";
						
			if($mut->{denovo} eq "yes") {
				
				# skip mutation if either parent is not sampled or not enough reference coverage
				next if ( (!exists $fatherCov{$key}) || (!exists $motherCov{$key}) );
				next if ( ($fatherCov{$key} < $min_cov) || ($motherCov{$key} < $min_cov) ); # min cov requirement
				
				if(!exists $denovoSVs{$k}) {
					$denovoSVs{$k} = $mut;
					$denovoL2K{$key} = $l2k->{$key};
				}
			}
			else { # inherithed or only in parents
				if($mut->{inheritance} ne "no") {
					#$vhash{$k} = $mut;
					if(!exists $inheritedSVs{$k}) {
						$inheritedSVs{$k} = $mut;
						$inheritedL2K{$key} = $l2k->{$key};
					}
				}
				else { # only in parents
					print STDERR "$key\n";
				}
			}
		}
	}
	
	#print denovos loc2keys table to file
	printLoc2key("$WORK/denovoloc2keys.txt", \%denovoL2K);
	printLoc2key("$WORK/inheritedloc2keys.txt", \%inheritedL2K);
	
	$denovo_dbm_obj->UnLock;
	$inherited_dbm_obj->UnLock;
}

## Compute best state for the location of the mutation
#####################################################
sub computeBestState {
	
	my $vars;	
	my $covthr = 0; # min cov threshold to remove sequencing errors
	
	foreach my $ID (@members) {
		
		if($ID eq "dad") { $vars = \%fatherSVs; }
		if($ID eq "mom") { $vars = \%motherSVs; }
		if($ID eq "sib") { $vars = \%sibSVs; }
		if($ID eq "aff") { $vars = \%selfSVs; }
		
		foreach my $currKey (keys %$vars) {
			#print STDERR "$currKey\n";
			
			next if($currKey eq "numPaths");
	
			my $var = $vars->{$currKey};
			my $chr = $var->{chr};
			my $pos = $var->{pos}; 
			my $loc = "$chr:$pos";
	
			my $refState;
			#my $mutState;
			my $altState;
			
			my $refCov;
			my $mutCov;
			my $altCov;
						
			foreach my $role (@members) {
				#print STDERR "$role: ";
				
				my $varh;	
				my $refhash;
				my $l2k;
				
				if($role eq "dad") { $varh = \%fatherSVs; $refhash = \%fatherCov; $l2k = \%fatherL2K; }
				if($role eq "mom") { $varh = \%motherSVs; $refhash = \%motherCov; $l2k = \%motherL2K; }
				if($role eq "sib") { $varh = \%sibSVs; $refhash = \%sibCov; $l2k = \%sibL2K; }
				if($role eq "aff") { $varh = \%selfSVs; $refhash = \%selfCov; $l2k = \%selfL2K; }
								
				$refState->{$role} = 0; $altState->{$role} = 0;
				$refCov->{$role} = 0; $mutCov->{$role} = 0; $altCov->{$role} = 0;
				
				my $sum = 0;
				my $cnt = 0;
				if ( exists($refhash->{$loc}) ) { # if there is reference coverage at location
					## compute reference coverage at locus
					for (my $t = $pos; $t<($pos+$kmer+1); $t++) { # compute avg coverage over window of size K
						my $loc2 = "$chr:$t";
						if ( exists $refhash->{$loc2} && 
							!( exists($l2k->{$loc2}) && ( (scalar @{$l2k->{$loc2}->{keys}})!=0) ) ) {
							$sum += $refhash->{$loc2};
							$cnt++;
						}					
					}
				}
				my $RCov = $sum;
				if($cnt>0) { $RCov = floor($sum/$cnt); }
				if($RCov > $covthr) { $refState->{$role} = 1; $refCov->{$role} = $RCov; }
				else { $refState->{$role} = 0; $refCov->{$role} = 0; }
		
				## compute coverage of mut at locus
				my $ACov = 0;
				if( exists($varh->{$currKey}) ) {
					my $altMut = $varh->{$currKey};
					$ACov = $altMut->{mincov};
				}
				if($ACov > $covthr) { $altState->{$role} = 1; $mutCov->{$role} = $ACov; }
				else { $altState->{$role} = 0; $mutCov->{$role} = 0; }
				
				## compute alternative allele coverage
				my $OCov = 0;
				if( exists($l2k->{$loc}) ) {
					foreach my $altKey (@{$l2k->{$loc}->{keys}}) {
						if( exists($varh->{$altKey}) && ($altKey ne $currKey) ) {
							my $altMut = $varh->{$altKey};
							$OCov += $altMut->{mincov};
						}
					}
				}
				if($OCov > $covthr) { $altCov->{$role} = $OCov; }
				else { $altCov->{$role} = 0; }
			}
			
			# 2 2 2 1/0 0 0 1 (ref: M,F,A,S / qry: M,F,A,S)
			my @RS=(); my @AS=();
			my @RC=(); my @AC=(); my @OC=();
			foreach my $role (@members) { # the order or members must be M,F,A,S
				if( ($refState->{$role} == 1) && ($altState->{$role} == 1) ) {
					push @RS , 1;
					push @AS , 1;
				}
				elsif( ($refState->{$role} == 1) && ($altState->{$role} == 0) ) {
					push @RS , 2;
					push @AS , 0;
				}
				elsif( ($refState->{$role} == 0) && ($altState->{$role} == 1) ) {
					push @RS , 0;
					push @AS , 2;
				}
				elsif( ($refState->{$role} == 0) && ($altState->{$role} == 0) ) {
					push @RS , 0;
					push @AS , 0;
				}
				
				push @RC , $refCov->{$role};
				push @AC , $mutCov->{$role};
				push @OC , $altCov->{$role};
			}
			
			my $stateStr = "";
			my $covStr = "";
			#if($famType eq "quad") {
				$stateStr = "$RS[0]$RS[1]$RS[2]$RS[3]/$AS[0]$AS[1]$AS[2]$AS[3]";
				$covStr   = "$RC[0] $RC[1] $RC[2] $RC[3]/$AC[0] $AC[1] $AC[2] $AC[3]/$OC[0] $OC[1] $OC[2] $OC[3]";
			#}
			#elsif($famType eq "trio") {
			#	$stateStr = "$R[0] $R[1] $R[2]/$A[0] $A[1] $A[2]";	
			#	$covStr   = "$RC[0] $RC[1] $RC[2]/$AC[0] $AC[1] $AC[2]/$OC[0] $OC[1] $OC[2]";		
			#}
			
			$var->{bestState} = $stateStr;
			$var->{covState} = $covStr;
			
			# compute hom/het state and chi2score
			parseCovState($var, $ID);
			
			$vars->{$currKey} = $var;
		}
	}
}

## parse covState and compute hom/het state and chi2score
##########################################

sub parseCovState {

	my $sv = $_[0];
	my $role = $_[1];
	my $BS = $sv->{covState};
	
	my ($REF,$ALT,$OTH) = split("/",$BS);
	my @R = split(" ", $REF); # reference
	my @A = split(" ", $ALT); # alternative
	my @O = split(" ", $OTH); # other

	my $id = 0;
	if($role eq "mom") { $id = 0; }
	if($role eq "dad") { $id = 1; }	
	if($role eq "aff") { $id = 2; }
	if($role eq "sib") { $id = 3; }
	
	my $alleleCnt = 0;	
	if($R[$id] > 0) { $alleleCnt++; }
	if($A[$id] > 0) { $alleleCnt++; }
	if($O[$id] > 0) { $alleleCnt++; }
	
	my $totCov = $R[$id] + $A[$id] + $O[$id];
	$sv->{altcov} = $totCov - $sv->{mincov};

	if($alleleCnt > 1) { $sv->{zygosity} = "het"; }
	else { $sv->{zygosity} = "hom"; }
	
	my $covRatio = 0;
	my $chi2Score = "na";
	if ($totCov > 0) { 
		$covRatio = sprintf('%.2f', ($sv->{mincov} / $totCov) );
		$sv->{covratio} = $covRatio; 

		#Chi-squared Test for allele coverage
		my $o1 = $sv->{mincov}; # observed 1
		my $o2 = $sv->{altcov}; # observed 2
		my $e1 = $totCov/2; # expected 1
		my $e2 = $totCov/2; # expected 2

		my $term1 = (($o1-$e1)*($o1-$e1))/$e1;
		my $term2 = (($o2-$e2)*($o2-$e2))/$e2;

		$chi2Score = sprintf('%.2f', $term1 + $term2);		
		#if($currMut->{zygosity} eq "hom") { $chi2Score = 0; }
	}
	$sv->{chi2score} = $chi2Score;	
}

## parse bestState and update mutation info
##########################################
sub parseBestStateQuad {

	my $sv = $_[0];
	my $BS = $sv->{bestState};
	
	$sv->{denovo} = "no";
	
	if( ($BS eq "1222/1000") || ($BS eq "0222/2000") ) { 
		#$sv->{note} = "only in mom"; 
		$sv->{inheritance} = "no";	
	}
	elsif( ($BS eq "2122/0100") || ($BS eq "2022/0200") ) { 
		#$sv->{note} = "only in dad"; 
		$sv->{inheritance} = "no";	
	}
	elsif( ($BS eq "1122/1100") || ($BS eq "0022/2200") ||
		   ($BS eq "1022/1200") || ($BS eq "0122/2100") ) { 
		#$sv->{note} = "only in parents"; 
		$sv->{inheritance} = "no";
	}	
	else {
		my ($REF,$ALT) = split("/",$BS);
		my @R = split("", $REF); # reference
		my @A = split("", $ALT); # alternative
		my $m=0; my $d=1; my $a=2; my $s=3;
		if ( ($A[$m]>0) && ($A[$d]>0) && # mom and dad have alternative allele
			 ( ($A[$a]>0) || ($A[$s]>0) ) ) { # either kid has alternative allele
			#$sv->{note} = "inherited from both parents";
			$sv->{inheritance} = "both";
		}
		elsif ( ($A[$m]>0) && ($A[$d]==0) && # only mom has alternative allele
			    ( ($A[$a]>0) || ($A[$s]>0) ) ) { # either kid has alternative allele
			#$sv->{note} = "inherited from mom";
			$sv->{inheritance} = "mom";
		}
		elsif ( ($A[$m]==0) && ($A[$d]>0) && # only dad has alternative allele
			    ( ($A[$a]>0) || ($A[$s]>0) ) ) { # either kid has alternative allele
			#$sv->{note} = "inherited from dad";
			$sv->{inh} = "dad";
		}
		elsif ( ($A[$m]==0) && ($A[$d]==0) &&  # mom and dad like reference
			    ( ($A[$a]>0) || ($A[$s]>0) ) ) { # either kid has alternative allele
			#$sv->{note} = "denovo";
			$sv->{denovo} = "yes";
			$sv->{inheritance} = "no";
		}
		else { 
			#$sv->{note} = "unknown state"; 
			$sv->{inheritance} = "no";	
		}
		
		# check which kid has the mutation
		if( ($A[$a]>0) && ($A[$s]>0) ) {  $sv->{id} = "both"; }
		elsif( ($A[$a]>0) && ($A[$s]==0) ) {  $sv->{id} = "aff"; }
		elsif( ($A[$a]==0) && ($A[$s]>0) ) {  $sv->{id} = "sib"; }
	}
}


## export family mutations to file
#####################################################

sub exportSVs {
	
	print STDERR "-- Exporting SVs to file\n";
	
	my $max_cov = 10000000;
	
	# export denovos 
	#runCmd("export denovos", "$exportTool --db $WORK/denovos.db --bed $BEDFILE --format annovar --type all --mincov $min_cov --maxcov $max_cov --covratio $outratio > $WORK/denovos.${min_cov}x.all.txt");
	#runCmd("export denovos", "$exportTool --db $WORK/denovos.db --bed $BEDFILE --format annovar --type snp --mincov $min_cov --maxcov $max_cov --covratio $outratio > $WORK/denovos.${min_cov}x.snp.txt");
	#runCmd("export denovos", "$exportTool --db $WORK/denovos.db --bed $BEDFILE --format annovar --type indel --mincov $min_cov --maxcov $max_cov --covratio $outratio > $WORK/denovos.${min_cov}x.indel.txt");

	my $command_denovo = "$exportTool ".
		"--db $WORK/denovos.db ".
		"--bed $BEDFILE ".
		"--format $outformat ".
		"--type indel ". 
		"--mincov $min_cov ". 
		"--maxcov $max_cov ".
		"--covratio $outratio";
	if($intarget) { $command_denovo .= " --intarget"; }
	if ($outformat eq "annovar") { $command_denovo .= " > $WORK/denovo.${min_cov}x.indel.annovar"; }
	elsif ($outformat eq "vcf") { $command_denovo .= " > $WORK/denovo.${min_cov}x.indel.vcf"; }
	
	print STDERR "Command: $command_denovo\n" if($VERBOSE);
	runCmd("export denovos", $command_denovo);
	
	# export inherited
	#runCmd("export inherited", "$exportTool --db $WORK/inherited.db --bed $BEDFILE --format annovar --type all --mincov $min_cov --maxcov $max_cov --covratio $outratio > $WORK/inherited.${min_cov}x.all.txt");
	#runCmd("export inherited", "$exportTool --db $WORK/inherited.db --bed $BEDFILE --format annovar --type snp --mincov $min_cov --maxcov $max_cov --covratio $outratio > $WORK/inherited.${min_cov}x.snp.txt");
	#runCmd("export inherited", "$exportTool --db $WORK/inherited.db --bed $BEDFILE --format annovar --type indel --mincov $min_cov --maxcov $max_cov --covratio $outratio > $WORK/inherited.${min_cov}x.indel.txt");

	my $command_inh = "$exportTool ".
		"--db $WORK/inherited.db ".
		"--bed $BEDFILE ".
		"--format $outformat ".
		"--type indel ". 
		"--mincov $min_cov ". 
		"--maxcov $max_cov ".
		"--covratio $outratio";
	if($intarget) { $command_inh .= " --intarget"; }
	if ($outformat eq "annovar") { $command_inh .= " > $WORK/inherited.${min_cov}x.indel.annovar"; }
	elsif ($outformat eq "vcf") { $command_inh .= " > $WORK/inherited.${min_cov}x.indel.vcf"; }
	
	print STDERR "Command: $command_inh\n" if($VERBOSE);
	runCmd("export inherited", $command_inh);
}

## process family
#####################################################

sub processFamily {	
	callSVs(); # call mutations on each family member
	loadHashes(); # load mutations and coverage info
	computeBestState(); # compute bestState for each mutation
	findDenovos(); # detect denovo mutations in the kids
	exportSVs(); # export mutations to file
}

## do the job
##########################################

init();
printParams("$WORK/parameters.txt");
processFamily();

##########################################

my $time_taken = time - $start_time;

#if($VERBOSE) {
	elapsedTime($time_taken, "FindDenovos");
#}