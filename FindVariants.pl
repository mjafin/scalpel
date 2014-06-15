#!/usr/bin/perl

###################################################################
# FindVariants.pl
#
# Tool for detecting variants in one single exome 
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
use Utils;
use SequenceIO;
use HashesIO;
use Parallel::ForkManager;
use List::Util qw[min max];
#use MathRandom::Random qw(:all);
use Digest::MD5 qw(md5_hex);
#use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Getopt::Long;
use File::Spec;

use MLDBM::Sync;						# this gets the default, SDBM_File
use MLDBM qw(DB_File Storable);			# use Storable for serializing
use MLDBM qw(MLDBM::Sync::SDBM_File);	# use extended SDBM_File, handles values > 1024 bytes
use Fcntl qw(:DEFAULT);					# import symbols O_CREAT & O_RDWR for use with DBMs

$|=1;

# required arguments
my $bamfile;
my $bedfile;
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
my $sample = $defaults->{sample};
my $selected = $defaults->{selected};
my $outformat = $defaults->{format};
my $intarget;

#microassembler parameters:
my $dfs_limit = $defaults->{pathlimit};
my $max_tip_len = 25;
my $maxKmer = 90;

my $help = 0;
my $VERBOSE = 0;
my $COV2FILE = 0;

# programs via absolute path
#my $samtools = "samtools";
my $microassembler = "$Bin/Microassembler/Microassembler";
my $exportTool = "$Bin/ExportVariants.pl";
my $bamtools = "$Bin/bamtools-2.3.0/bin/bamtools";

my $rgfile = "readgroups.txt";

#my %pathways;
#my %exonshash;
my %variants;
my %exons;
my %locations;
my %readgroups;
my %refcov;
my %loc2keys;
###my %genome;
my $variants_dbm_obj;

my $argcnt = scalar(@ARGV);
my $start_time = time;

#####################################################
#
# Command line options processing
#
GetOptions(
	'help!'    => \$help,
	'verbose!' => \$VERBOSE,
	'cov2file' => \$COV2FILE,
    
	# required parameters
    'bam=s' => \$bamfile,
    'bed=s' => \$bedfile,
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
    'pathlimit=i'  => \$dfs_limit,
    'dir=s'        => \$WORK,
    'numprocs=i'   => \$MAX_PROCESSES,
    'sample=s'     => \$sample,
    'coords=s'     => \$selected,
    'format=s'     => \$outformat,

	# ouptut parameters
	'intarget!'    => \$intarget,
	'outratio=f'   => \$outratio,
    

) or usageSingle("FindVariants.pl");

#####################################################
#
# Initilization 
#
sub init()
{
	usageSingle("FindVariants.pl") if ($argcnt < 1);
	usageSingle("FindVariants.pl") if( $help );
	
	if( (!defined $bamfile) || (!defined $bedfile) || (!defined $REF) ) { 
		print STDERR "Required parameter missing!\n";
		usageSingle("FindVariants.pl");
	}
	
	mkdir $WORK if ! -r $WORK;
	
	if($MAX_PROCESSES == 1) { $MAX_PROCESSES = 0; }
}

#####################################################

sub printParams { 

	my $file = $_[0]; 
	
	print STDERR "-- Print parameters to $file\n";
	
	open PFILE, "> $file" or die "Can't open $file ($!)\n";
	
	print PFILE "Parameters: \n";
	print PFILE "-- k-mer size (bp): $kmer\n";
	print PFILE "-- coverage threshold: $cov_threshold\n";
	print PFILE "-- low coverage threshold: $tip_cov_threshold\n";
	print PFILE "-- size (bp) of the left and right extensions (radius): $radius\n";
	print PFILE "-- window-size size (bp): $windowSize\n";
	print PFILE "-- step size (bp): $delta\n";
	print PFILE "-- minimum mapping quality: $map_qual\n";
	print PFILE "-- minimum coverage for mutation: $min_cov\n";
	print PFILE "-- minimum coverage ratio for sequencing errors: $covratio\n";
	print PFILE "-- minimum coverage ratio for exporting mutation: $outratio\n";
	print PFILE "-- max number of mismatches for near-perfect repeats: $maxmismatch\n";
	print PFILE "-- limit number of sequence paths: $dfs_limit\n";
	print PFILE "-- output directory: $WORK\n";
	print PFILE "-- max number of parallel jobs: $MAX_PROCESSES\n";
	print PFILE "-- bam file: $bamfile\n";
	print PFILE "-- bed file: $bedfile\n";
	print PFILE "-- reference file: $REF\n";
	print PFILE "-- sample: $sample\n";
	print PFILE "-- output format for variants: $outformat\n";
	print PFILE "-- file of selected coordinates: $selected\n";
	if($intarget) { print PFILE "-- output variants in target? yes\n"; }
	else { print PFILE "-- output variants in target? no\n"; }
	
	close PFILE;
}

## build directory structure and run assembly
#####################################################

sub run {

	my $fdir = "$WORK";
	
	if (! -r $fdir) { mkdir $fdir; }

	printParams("$fdir/parameters.txt");
		
	if($selected ne "null") {
		loadCoordinates("$selected", \%exons, $VERBOSE);
	}
	else {
		loadExonsBed("$bedfile", \%exons, $radius, $VERBOSE);
	}
	
	#load genome from fasta file
	###loadGenomeFasta($REF, \%genome);
	
	if (-e "$fdir/variants.db.dir") { runCmd("remove database files", "rm $fdir/variants.db.*"); }	
	$variants_dbm_obj = tie %variants, 'MLDBM::Sync', "$fdir/variants.db", O_CREAT|O_RDWR, 0640;	
	#$variants_dbm_obj->SyncCacheSize('1000M');    # 1000 Megabyte max memory used
	
	#print "$bamfile\n";
	
	# get absolute path
	my $bam_abs_path = File::Spec->rel2abs($bamfile);

	## link the input files	
    if (-e $bam_abs_path) { 
		symlink "$bam_abs_path", "$fdir/bamfile.bam"; 
	}
    else { 
		print STDERR "Can't find bamfile: $bamfile\n";
		next;
    }

    if (-e "$bam_abs_path.bai") {
		symlink "$bam_abs_path.bai", "$fdir/bamfile.bam.bai";
    }
    else {
		print STDERR "Indexing BAM file...\n";
		runCmd("index bamfile", "$bamtools index -in $bam_abs_path");
		symlink "$bam_abs_path.bai", "$fdir/bamfile.bam.bai";
    }

	parseHeader($WORK, \%readgroups);
	saveReadGroup("$sample", $WORK, \%readgroups, $rgfile);
	
	print STDERR "Assembly Exons\n";
	assemblyExons(\%variants, $variants_dbm_obj, \%refcov, \%loc2keys);

	if ($COV2FILE) { printRefCov("$WORK/refcov.txt", \%refcov); }
	printLoc2key("$WORK/loc2keys.txt", \%loc2keys);
	
	# export variants to file
	exportSVs();
}

## export family mutations to file
#####################################################

sub exportSVs {

	#$exportvars --db $DIR/variants.db --bed $BEDFILE --format annovar --type all --mincov $mincov --maxcov $maxcov > $DIR/variants.txt	
	print STDERR "-- Exporting SVs to file\n";
	
	my $max_cov = 10000000;
	
	# export denovos 
	#runCmd("export denovos", "$exportTool --db $WORK/variants.db --bed $bedfile --format annovar --type all --mincov $min_cov --maxcov $max_cov --covratio $outratio > $WORK/variants.${min_cov}x.all.txt");
	#runCmd("export denovos", "$exportTool --db $WORK/variants.db --bed $bedfile --format annovar --type snp --mincov $min_cov --maxcov $max_cov --covratio $outratio > $WORK/variants.${min_cov}x.snp.txt");
	#runCmd("export denovos", "$exportTool --db $WORK/variants.db --bed $bedfile --format annovar --type indel --mincov $min_cov --maxcov $max_cov --covratio $outratio > $WORK/variants.${min_cov}x.indel.txt");
	
	my $command = "$exportTool ".
		"--db $WORK/variants.db ".
		"--bed $bedfile ".
		"--format $outformat ".
		"--type indel ". 
		"--mincov $min_cov ". 
		"--maxcov $max_cov ".
		"--covratio $outratio";
	if($intarget) { $command .= " --intarget"; }
	if ($outformat eq "annovar") { $command .= " > $WORK/variants.${min_cov}x.indel.annovar"; }
	elsif ($outformat eq "vcf") { $command .= " > $WORK/variants.${min_cov}x.indel.vcf"; }
	
	print STDERR "Command: $command\n" if($VERBOSE);
	runCmd("ExportVariants", $command);
}

#####################################################


## Write the list of read groups per family member on file 
## Extract the sequence around indel from reference
## Create input files for microassembler
## Run microassembler

sub assemblyExons {
	
	my $hash = $_[0];
	my $dbm_obj = $_[1];
	my $ref_hash = $_[2];
	my $loc2keys_hash = $_[3];
	
	my $count_exons = 0;
	my $count_regions = 0;
	my @regions;
	my @region_files;
	
	my $PREFIX = "$WORK";
						
	#foreach my $exonkey (sort {$a <=> $b} keys %exons) { # for each chromosome
	foreach my $exonkey (keys %exons) { # for each chromosome

		$count_exons = 0;
		#foreach my $exon (sort {$a->{start} <=> $b->{start}} @{$exons{$exonkey}}) { # for each exon
		foreach my $exon (@{$exons{$exonkey}}) { # for each exon
			
			$count_exons++;
			#last if($count_exons > 5);
		
			#my %pathways;
			#my %exonshash;
			
			# Forks and returns the pid for the child:
			#my $pid = $pm->start and next;
	
			my $start = $exon->{start};
			my $end = $exon->{end};
			my $chr = $exon->{chr};
			#print STDERR "$chr:$start-$end\n";
			
			## extend region left and right by radius
			my $left = $start-$radius;
			if ($left < 0) { $left = $start; }
			my $right = $end+$radius;
			#if ($right > $end) { $right = $end; }
			#print STDERR "$chr:$left-$right\n";
			
			my $exonSize = $right-$left;
							
			# generate list of regions to examine			
			my $l = $left;
			my $u = $l+$windowSize;
			if($u >= $right) { 
				$u = $right;
				my $region;	
				$region->{chr} = $chr;				
				$region->{left} = $l;
				$region->{right} = $u;
				push @regions, $region;
				#print "$region->{chr}:$region->{left}-$region->{right}\n";
				$count_regions++;
			}
			else {
				while($u < $right) {
				
					my $region;
					$region->{chr} = $chr;		
					$region->{left} = $l;
					$region->{right} = $u;
					$region->{exonstart} = $left;
					$region->{exonend} = $right;					
					push @regions, $region;
					$count_regions++;
			
					$l += $delta;
					$u += $delta;
					if($u > $right) { 
						$u = $right; 
						$l = $u-$windowSize; 
						
						my $region;
						$region->{chr} = $chr;
						$region->{left} = $l;
						$region->{right} = $u;
						$region->{exonstart} = $left;
						$region->{exonend} = $right;
						
						push @regions, $region;
						#print "$region->{chr}:$region->{left}-$region->{right}\n";
						$count_regions++;
					}
				}
			}
		}
	}
	
	#update statistics in the database
	if(exists $hash->{stats})  { 
		my $sts = $hash->{stats}; 
		$sts->{num_exceptions} = 0;
		$sts->{num_dfs_limit} = 0;
		$sts->{num_partial_align} = 0;
		$sts->{num_with_cycles} = 0;
		$sts->{num_ok} = 0;
		$hash->{stats} = $sts;
	}
	
	#schedule jobs
	schedule_jobs(\@regions, $hash, $dbm_obj, $ref_hash, $loc2keys_hash, "regions");

	my $stats = $hash->{stats}; 

	#if ($VERBOSE) {
	if(exists $hash->{stats})  {
		my $sts = $hash->{stats};
		print STDERR "\nDfs limit: $sts->{num_dfs_limit}\n";
		print STDERR "Partial alignments: $sts->{num_partial_align}\n";
		print STDERR "Cycles: $sts->{num_with_cycles}\n";
		print STDERR "OK assemblies: $sts->{num_ok}\n";
		print STDERR "Exceptions: $sts->{num_exceptions}\n";
		print STDERR "Repeats: $sts->{num_repeats}\n\n";
	}
	#}
}

# partition the region to assemble into equal size 
# random sets and run the parallel jobs
#####################################################

sub schedule_jobs {
	
	my $array = $_[0];
	my $hash = $_[1];
	my $dbm_obj = $_[2];
	my $ref_hash = $_[3];
	my $loc2keys_hash = $_[4];
	my $fprefix = $_[5];
	
	my $PREFIX = "$WORK";
	
	# random permutation of the regions to assemble
	#my @regions = random_permutation(@$array);
	my @regions = @$array;
	my $arraySize = $#regions + 1;

	print STDERR "start assembly of $arraySize regions.\n"; 

	# partition regions
	my $step = $arraySize;
	if ($MAX_PROCESSES>0) {
		$step = ceil($arraySize/$MAX_PROCESSES);
	}	
	print STDERR "stepSize: $step\n";

	for(my $i = 0; $i < $arraySize; $i+=$step) {
	
		my $s = $i;
		my $e = min( ($arraySize - 1), ($i + $step-1) );
		my $ref_file = "$fprefix-$s-$e.fa";
		#push @region_files, $ref_file;
			
		open FILE, "> $PREFIX/$ref_file" or die "Can't open $PREFIX/$ref_file ($!)\n";
		foreach my $reg (@regions[$s..$e]) { # for each region
			my $chr = $reg->{chr};
			my $left = $reg->{left};
			my $right = $reg->{right};
			my $exonstart = $reg->{exonstart};
			my $exonend = $reg->{exonend};
			
			#my $chrom = $chr;
			#if ($chr =~ /^chr/) { $chrom = substr($chr,3); }

			print FILE ">$chr:$left-$right\n";
			### die "Undefined sequence ($chr)\n" if (!exists($genome{$chr}));
			my ($header, $seq) = split(/\n/, `samtools faidx $REF $chr:$left-$right`, 2);
			$seq =~ s/[\n\r\s]+//g; 
			###my $seq = substr($genome{$chr}->{seq}, $left-1, $right-$left+1);
			#my $seq = substr($genome{$chr}->{seq},$exonstart, $exonend-$exonstart+1);
			for(my $i=0; $i<length($seq); $i+=80) {
				my $str = substr($seq,$i,80);
				print FILE "$str\n";
			}
		}
		close FILE;
	}
	
	# free genome memory before parallelization
	###for my $k (keys %genome) { delete $genome{$k}; }
	
	my $pm = new Parallel::ForkManager($MAX_PROCESSES);
	my $count = 0;
	for(my $i = 0; $i < $arraySize; $i+=$step) {
		$count++;
	
		my $pid = $pm->start and next;

		my $s = $i;
		my $e = min( ($arraySize - 1), ($i + $step-1) );
		my $ref_file = "$fprefix-$s-$e.fa";
	
		print STDERR "$count. [$s..$e]\n";
		
		assemblyRegion($count, $ref_file, \%$hash, $dbm_obj);
		$pm->finish; # Terminates the child process
	}
	$pm->wait_all_children; # wait for all parallel jobs to complete

	# load reference coverage info from files
	loadRefCovMerge($ref_hash, $count, $WORK);
	loadLoc2KeyMerge($loc2keys_hash, $ref_hash, $count, $kmer, $WORK);
	
	# compute zygosity for all mutations
	computeHomHetState($hash, $ref_hash, $loc2keys_hash);
}

## Compute homozygous vs heterozygous state and  
## flags possible sequencing errors
#####################################################

sub computeHomHetState {
	
	print STDERR "Computing hom/het state...\n";
	
	my $vars = $_[0];
	my $refh = $_[1];
	my $l2k  = $_[2];
	
	my $hom_cnt = 0;
	my $het_cnt = 0;
	
	# improve performance: tie once to database and do all read/wrote operations.
	$variants_dbm_obj->Lock;

	foreach my $currKey (keys %$vars) {
		
		next if($currKey eq "stats"); # skip stats info
				
		my $currMut = $vars->{$currKey};
		my $chr = $currMut->{chr};
		my $pos = $currMut->{pos};
		
		my $loc = "$chr:$pos";
		
		my $alleleCnt = 0;
		my $Rcov = 0;
		my $Acov = 0;
		my $Ocov = 0;
		my $totCov = 0;
		
		# reference coverage
		my $sum = 0;
		my $cnt = 0;
		if ( (exists $refh->{$loc}) ) { # if cov at location					
			## compute reference coverage at locus
			for (my $t = $pos; $t<($pos+$kmer+1); $t++) { # compute avg coverage over window of size K
				my $loc2 = "$chr:$t";
				if ( (exists $refh->{$loc2})  && 
					!( (exists $l2k->{$loc2}) && ( (scalar @{$l2k->{$loc2}->{keys}})!=0 ) ) ) {
					$sum += $refh->{$loc2};
					$cnt++;
				}					
			}
		}
		$Rcov=$sum;
		if($cnt>0) { $Rcov = floor($sum/$cnt); }
		
		# alternative allele coverage
		$Acov = $currMut->{mincov};
		
		# other coverage
		$Ocov = 0;
		if(exists($l2k->{$loc})) {
			foreach my $altKey (@{$l2k->{$loc}->{keys}}) {
				if( exists($vars->{$altKey}) && ($altKey ne $currKey) ) {
					my $altMut = $vars->{$altKey};
					$Ocov += $altMut->{mincov};
				}
			}
		}

		$totCov = $Rcov + $Acov + $Ocov;		
		$currMut->{altcov} = $totCov - $currMut->{mincov};
		
		if($Rcov > 0) { $alleleCnt++; }
		if($Acov > 0) { $alleleCnt++; }
		if($Ocov > 0) { $alleleCnt++; }
		
		if($alleleCnt > 1) {
			$currMut->{zygosity} = "het";
			$het_cnt++;
		}
		else { # alleCnt <= 1
			$currMut->{zygosity} = "hom";
			$hom_cnt++;
		}
		
		# update mutation in db
		$vars->{$currKey} = $currMut;
	}
	
	$variants_dbm_obj->UnLock;

	print STDERR "Num. Homozygous: $hom_cnt\n";
	print STDERR "Num. Heterozygous: $het_cnt\n";
}

# assembly a specific region
#####################################################

sub assemblyRegion {

	my $num = $_[0];	
	my $refs_file = $_[1];
	my $hash = $_[2];
	my $dbm_obj = $_[3];
	
	#my %pathways;
	#my %exonshash;
	
	my $PREFIX = "$WORK";
	my $log_dir = "$WORK/logs";
	if (! -r $log_dir) { mkdir $log_dir; }
	
	my $outfile = "$log_dir/$refs_file.log";
	
	my $command = "$microassembler -R -I -C ";
	if($VERBOSE) { $command .= "-v "; }
	$command .= "-F $dfs_limit ".
		"-k $kmer ".
		"-K $maxKmer ". 
		"-M $maxmismatch ".
		"-l $max_tip_len ". 
		"-c $cov_threshold ".
		"-x $covratio ".
		"-d $tip_cov_threshold ". 
		"-b $map_qual ". 
		"-p $PREFIX ".
		"-m bamfile.bam ".
		"-g $rgfile ".
		" -r $refs_file ".
		">& $outfile";
		
	## run microassembler on each fasta file 
	runCmd("microassembly", "$command"); 
	runCmd("gzip", "gzip -f $outfile"); 
	
	parseLogFile("$outfile", $num, \%$hash, $dbm_obj);
			
	# remove reference file
	runCmd("delete file", "rm $PREFIX/$refs_file");
}

## parse the log file generated by microassembler and populate data structures containing the output statistics
#####################################################

sub parseLogFile {

	my $logFile = $_[0];
	my $num = $_[1];
	my $hash = $_[2];
	my $dbm_obj = $_[3];

	#my $PATHS = new IO::Uncompress::Gunzip "$logFile.gz" or die "IO::Uncompress::Gunzip failed: $GunzipError\n";
	open PATHS, "gunzip -c $logFile.gz|" or die "Can't open $logFile.gz ($!)\n!";

	my $curc;
	my $curs;
	my $cure;

	my @pathcov;
	my @pathdiff;

	my %tmp_hash;
	my %tmp_ref_cov;
	my %tmp_loc2keys;
	
	my $stats;
	$stats->{num_dfs_limit} = 0;
	$stats->{num_partial_align} = 0;
	$stats->{num_with_cycles} = 0;
	$stats->{num_ok} = 0;
	$stats->{num_exceptions} = 0;
	$stats->{num_repeats} = 0;

	my $exon;
	my @pathways;

	while (<PATHS>)
	{
		if (/== Processing\s+\d+:\s+([\w\.]+):(\d+)-(\d+)/) {
			
			$curc = $1; # chr
			$curs = $2; # start
			$cure = $3; # end
		}
		elsif (/cov: (\S+)/) {
			my $cov = $1;

			$exon->{start}  = $curs;
			$exon->{end}    = $cure;
			$exon->{cov}    = $cov;
			$exon->{status} = "-";

			my $exonkey = sprintf("%s:%d:%d", $curc, $curs, $cure);

			#$exonshash->{$exonkey} = $exon;
		}
		elsif (/ref trim5: (\d+) trim3: (\d+) uncovered: (\d+) ref_dist: (\d+)/) {
			$exon->{trim5} = $1;
			$exon->{trim3} = $2;
			$exon->{uncovered} = $3;
			$exon->{ref_dist} = $4;
		}
		elsif (/allcycles: (\d+)/) {
			$exon->{cycles} = $1;
			if ($exon->{cycles} >= 1) { $stats->{num_with_cycles} += 1; }			
		}
		elsif (/exception/) {
			$exon->{status} = "excp";
			$stats->{num_exceptions} += 1;
		}
		elsif (/Found repeat in/ || /Found cycle in/) {
			$exon->{status} = "repeat";
			$stats->{num_repeats} += 1;
		}
		elsif (/WARNING: DFS_LIMIT/) {
			$exon->{status} = "dfs";
			$stats->{num_dfs_limit} += 1;
		}
		elsif(/^>sp/) {
			$exon->{status} = "partial";
			$stats->{num_partial_align} += 1;
		}
		elsif(/^c':/) {
			chomp;
			my $covstr = (split(/:/, $_))[1];
			@pathcov = split(/ /, $covstr);
		}
		elsif(/^d':/) {
			chomp;
			my $diffstr = (split(/:/, $_))[1];
			@pathdiff = split(//,$diffstr);
		}
		elsif (/^>p_([\w\.]+):(\d+)-(\d+)_(\d+)/) {
			$exon->{status} = "ok";
			my $chr = $1;

			my $path;
			$path->{chr} = $chr;
			$path->{rstart} = $2;
			$path->{rend}   = $3;
			$path->{trim5}  = $exon->{trim5};
			$path->{trim3}  = $exon->{trim3};	
			$path->{pcnt}   = $4;
			@{$path->{cov}} = @pathcov;
			@{$path->{diff}} = @pathdiff;
			chomp;

			# >p_9:96439830-96440030_1 cycle: 0 match: 200 snp: 0 ins: 0 del: 1 96439931:T|-|25.7

			my ($tag, $cc, $cycle, $mm, $match, $ss, $snp, $ii, $ins, $dd, $del, $info) = split(/\s/, $_, 12);

			$path->{cycle} = $cycle;
			$path->{match} = $match;
			$path->{snp} = $snp;
			$path->{ins} = $ins;
			$path->{del} = $del;
			$path->{info} = $info;

			#my $pathkey = sprintf("%s:%d:%d", $chr, $2, $3);
			if($path->{cycle} == 0) { ## skip paths with cycles
				push @pathways, $path;
			}
		}
		elsif (/FINISHED/) { # end of current region
			if((scalar @pathways) > 0) { 
				if( ($exon->{status} eq "ok") && ($exon->{cycles} == 0) ) { # skip regions hard to assemble and paths with at least one cycle
					# at this point the exon status is ok
					$stats->{num_ok} += 1;
					extractVariantsFromPath(\@pathways, $exon, \%tmp_hash, \%tmp_ref_cov, $stats, \%tmp_loc2keys);
				}
				undef (@pathways); # free up the memory resources occupied by the array
			}
		}
	}
	
	#copy mutations to database
	$dbm_obj->Lock;
	foreach my $key (keys %tmp_hash) {
		my $mut = $tmp_hash{$key};
		
		if( !(exists $hash->{$key}) ) {
			$hash->{$key} = $mut;
		}
		else { # denovo already in the table
			# update the coverage to the highest value found so far
			my $mut_old = $hash->{$key};
			if ($mut_old->{avgcov} < $mut->{avgcov}) { $mut_old->{avgcov} = $mut->{avgcov}; } # max avg coverage
			if ($mut_old->{mincov} < $mut->{mincov}) { $mut_old->{mincov} = $mut->{mincov}; } # max min coverage
			$hash->{$key} = $mut_old;
		}
	}
	
	#update statistics in the database
	if( !(exists $hash->{stats}) ) { 
		$hash->{stats} = $stats; 
	}
	else { 
		my $old_stats = $hash->{stats}; 
		$old_stats->{num_repeats} += $stats->{num_repeats};
		$old_stats->{num_exceptions} += $stats->{num_exceptions};
		$old_stats->{num_dfs_limit} += $stats->{num_dfs_limit};
		$old_stats->{num_partial_align} += $stats->{num_partial_align};
		$old_stats->{num_with_cycles} += $stats->{num_with_cycles};
		$old_stats->{num_ok} += $stats->{num_ok};
		$hash->{stats} = $old_stats;
	}
	$dbm_obj->UnLock;
	
	# write reference coverage to file	
	my $outfile = "$WORK/refcov.$num.txt";
	open FILE, "> $outfile" or die "Can't open $outfile ($!)\n";
	foreach my $key (keys %tmp_ref_cov) {
		print FILE "$key\t$tmp_ref_cov{$key}\n";
	}
	close FILE;
	
	# write mutations by positions to file	
	$outfile = "$WORK/loc2keys.$num.txt";
	open FILE, "> $outfile" or die "Can't open $outfile ($!)\n";
	foreach my $loc (keys %tmp_loc2keys) {
		print FILE "$loc";		
		foreach my $key (@{$tmp_loc2keys{$loc}}) { print FILE "\t$key"; }
		print FILE "\n";
	}
	close FILE;

	close PATHS;
}

## extract mutations from path
##########################################

sub extractVariantsFromPath {
	
    my $pathways = $_[0];
    my $exon = $_[1];
	my $hash = $_[2];
	my $ref_cov_hash = $_[3];
	my $stats = $_[4];
	my $mutbypos_hash = $_[5];

	my $outflag = 0;
	
	my $exoncov = $exon->{cov};
	my $exonsta = $exon->{status};
	my $exoncycle = $exon->{cycles};
	$exoncycle = 0 if !defined $exoncycle;

	foreach my $path (@{$pathways}) { ## for each path in the region
		
		updateRefCoverage($path, $ref_cov_hash); # updated coverage info for matching reference

		my $chr = $path->{chr};	
		#if ($chr =~ /^chr/) { $chr = substr($chr,3); }
		my $info = $path->{info};
		$info = "" if !defined $info;
		#print "---$info---\n";
		my @var = split /\s+/, $info;
					
		foreach my $v (@var) { ## for each variant
			my ($pos, $ref, $qry, $avgcov, $mincov, $prevbpref, $prevbpalt) = split /[:|]/, $v;
			
			my $t = "";
			my $l = 0;
			if   ( ($ref !~ /-+/) && ($qry =~ /-+/) ) {$t = "del"; $l = length($qry); } # deletion
			elsif( ($ref =~ /-+/) && ($qry !~ /-+/) ) {$t = "ins"; $l = length($qry); } # insertion
			elsif( ($ref !~ /-+/) && ($qry !~ /-+/) ) {$t = "snp"; $l = length($qry); } # SNP

			$outflag = 1;
					
			my $mut;
			$mut->{chr}  = $chr;							
			$mut->{pos}  = $pos;
			$mut->{type} = $t;
			$mut->{len}  = $l;
			$mut->{ref}  = $ref;
			$mut->{seq}  = $qry;
			$mut->{prevbpref} = $prevbpref;
			$mut->{prevbpalt} = $prevbpalt;
			$mut->{avgcov} = $avgcov;
			$mut->{mincov} = $mincov;
			$mut->{status} = $exonsta;
			$mut->{zygosity} = "na";
			$mut->{altcov} = 0;
			$mut->{inheritance} = "na";
		
			my $ref_encoded = md5_hex($ref);
			my $qry_encoded = md5_hex($qry);
					
			my $key = sprintf("%s:%d:%s:%d:%s:%s", $chr, $pos, $t, $l, $ref_encoded, $qry_encoded);
			my $loc = sprintf("%s:%d", $chr, $pos);
			#print "$key\n";

			## add mutation to the hash table
			if( !(exists $hash->{$key}) ) {
				$hash->{$key} = $mut;
				push @{$mutbypos_hash->{$loc}}, $key;
			}
			else { # already in the table
				# update the coverage to the highest value found so far
				my $mut_old = $hash->{$key};
				if ($mut_old->{avgcov} < $mut->{avgcov}) { 
					$mut_old->{avgcov} = $mut->{avgcov}; # max avg coverage
				} 
				if ($mut_old->{mincov} < $mut->{mincov}) { 
					$mut_old->{mincov} = $mut->{mincov}; # max min coverage
				} 
				$hash->{$key} = $mut_old;
			}
		}
	}
	#return 1 if at least one event is detected
	return $outflag;
}

## update covearge info for bases matching the reference
##########################################

sub updateRefCoverage {
	
	my $p = $_[0];
	my $ref_hash = $_[1];
	
	my @cov = @{$p->{cov}};
	my @diff = @{$p->{diff}};
  	my $chr = $p->{chr};
  	my $start = $p->{rstart};
  	my $trim5 = $p->{trim5};
  	#my $end = $p->{rend};
	#my $L = $end-$start+1; # length of the region
	my $len = scalar @diff;
		
	my $r = $start+$trim5; # reference index
	my $j = 0; # used to scan the coverage array
	for (my $i = 0; $i<$len; $i++) {

		my $key = "$chr:$r";

		if ($diff[$i] eq " ") { # if matching the reference
			#if( ($i>0) && ($diff[$i-1] eq " ") ) { # test needed to correctly handle homozygous ins
			if( ($i>0) && ($diff[$i-1] ne "^") ) { # test needed to correctly handle homozygous ins
				if ( (exists $ref_hash->{$key}) ) { # already in the table
					my $cov = $ref_hash->{$key};
					if ($cov < $cov[$j]) { # keep max coverage
						$cov = $cov[$j]; 
					}
					$ref_hash->{$key} = $cov;
				}
				else { # not in the table
					$ref_hash->{$key} = $cov[$j];
				}
			}
			$j++; $r++;
		}
		elsif ($diff[$i] eq "x") { $j++; $r++; } # snp 
		elsif ($diff[$i] eq "^") { $j++; } # ins 
		elsif ($diff[$i] eq "v") { $r++; } # del
	}	
}

## do the job
##########################################

init();
run();

##########################################

my $time_taken = time - $start_time;

#if($VERBOSE) {
	elapsedTime($time_taken, "FindVariants");
#}
