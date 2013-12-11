package SequenceIO;

###################################################################
# SequenceIO
#
# Package for common IO routines
# 
#  Author: Giuseppe Narzisi 
#    Date: December 11, 2013
#
###################################################################

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(loadCoordinates loadExonsBed loadGenomeFasta parseHeader saveReadGroup);

use strict;
use warnings;
use FindBin qw($Bin);
use lib $Bin; # add $Bin directory to @INC
use Utils;

my $bamtools = "$Bin/bamtools-2.3.0/bin/bamtools";

## load selected location (coordinates)
#####################################################
sub loadCoordinates {
	
	my $file = $_[0];
	my $exons = $_[1];
	my $VERBOSE = $_[2];
	
	print STDERR "Loading coordinates...";
	
	if ($file ne "null") {
		 
		open SELECTED, "< $file" or die "Can't open $file ($!)\n";

		my $cntall_locs = 0;
		while (<SELECTED>) {
			chomp;
			next if ($_ =~ /^#/); # skip comments
			
			my ($chr, $pos) = split /\t/, $_, 2;
			#if ($chr =~ /^chr/) { $chr = substr($chr,3); }
			
			my $exon;
			$exon->{chr} = $chr;
			$exon->{start} = $pos;
			$exon->{end} = $pos;
			
			push @{$exons->{$chr}}, $exon;
			$cntall_locs++;
		}
		close SELECTED;

		#if($VERBOSE) {
			print STDERR "$cntall_locs locations.\n";
		#}
	}
}

## load exons list
#####################################################
sub loadExonsBed {
	
	print STDERR "Loading targets from BED file...";

	my $file = $_[0];
	my $exons = $_[1];
	my $radius = $_[2];
	my $VERBOSE = $_[3];
	
	open EXONSLIST, "< $file" or die "Can't open $file ($!)\n";
	
	my $cntall_exons = 0;
	my $cntovl_exons = 0;

	my $prev_exon;
	$prev_exon->{chr} = 0;
	$prev_exon->{start} = 0;
	$prev_exon->{end} = 0;
	while (<EXONSLIST>)
	{
		chomp;
		next if ($_ =~ /^#/); # skip comments
		
		#my ($chr, $start, $end) = split /\t/, $_, 3;
		my @array = split /\t/, $_;
		my $chr = $array[0];
		my $start = $array[1];
		my $end = $array[2];
		#if ($chr =~ /^chr/) { $chr = substr($chr,3); }

		# exons coordinate are sorted in the bed file
		if ( ($prev_exon->{chr} eq $chr) && ($start <= $prev_exon->{end}) ) { 
			#print STDERR "exons regions overlap: [$prev_exon->{start},$prev_exon->{end}] [$start,$end]\n";
			$prev_exon->{end} = $end;
			$cntovl_exons++;
		}
		else {
			my $exon;
			$exon->{chr} = $chr;
			$exon->{start} = $start;
			$exon->{end} = $end;
			
			## extend region left and right by radius
			#my $l = $start-$radius;
			#my $u = $end+$radius;
	
			push @{$exons->{$exon->{chr}}}, $exon;
			
			# update prev exon
			$prev_exon = $exon;
			
			$cntall_exons++;
		}
		#last if($cntall_exons > 10);
	}
	
	if($VERBOSE) {
		# print exons
		foreach my $k (keys %$exons) { # for each chromosome
			foreach my $exon (@{$exons->{$k}}) { # for each exon
				print STDERR "$exon->{chr}\t$exon->{start}\t$exon->{end}\n";
			}
		}
	}

	close EXONSLIST;

	#if($VERBOSE) {
		print STDERR "$cntall_exons targets (filtered $cntovl_exons overlapping).\n";
	#}
}

# load the genome in fasta format
##########################################
sub loadGenomeFasta {
    
    print STDERR "Loading genome from FASTA file...";
    
    my $file = $_[0];
    my $genome = $_[1];
    
    open FASTAFILE, "< $file" or die "Can't open $file ($!)\n";

    my $header = "";
    my $SDNA = "";
    my $chr = "";
	my $cnt = 0;
    while ( <FASTAFILE> ) {
        chomp;
        #print "$_\n";
        if($_ =~ m/>/) {
			$cnt++;
            $header = $_;
            #print "$header\n";
            if($chr ne "") {
                #print "$chr\n";
                $genome->{$chr}->{seq} = $SDNA;
            }
            
            $header =~ s/^>//; # remove ">"
            $header =~ s/\s+$//; # remove trailing whitespace

            my ($label, $tmp) = split / /, $header, 2;
			$chr = $label;			
			#if ($label =~ /^chr/) { $chr = substr($label,3); } # update chromosome label
            $SDNA = "";# clear out old sequence
        }
        else { 
            s/\s+//g; # remove whitespace
            $SDNA .= $_; 
        }
    }
    close FASTAFILE;

    if($chr ne "") { # handle last sequence
        #print "$chr\n";
        $genome->{$chr}->{seq} = $SDNA;
    }  

	#if($VERBOSE) {
		print STDERR "$cnt sequences.\n";
	#}
}

## process SAM header to extract read groups info
#####################################################
sub parseHeader {
	 
	my $PREFIX = $_[0];
	my $readgroups = $_[1];
	my $headerfile = "header.txt";
	
	# erease readgroup information from previous family
	for (keys %$readgroups) { delete $readgroups->{$_}; }
	
	## extract header from SAM file
	runCmd("extract SAM header", "$bamtools header -in $PREFIX/bamfile.bam > $PREFIX/$headerfile");
	
	## read the list of read groups per sample from header
	open HEADER, "< $PREFIX/$headerfile" or die "Can't open $PREFIX/$headerfile ($!)\n";
	
	while (<HEADER>) {
		chomp;
		next if($_ eq ""); # skip empty string
		next if($_ !~ /^\@/); # skip non-header line
		
		my @records = split /\t/, $_;
	  	if ($records[0] eq "\@RG") {
			#@RG     ID:READGROUP    SM:SAMPLE
			
	  		#print join(' ', @records), "\n";
	  		
	  		my $sm;
	  		my $id;
	  		foreach my $tag (@records) {
	  			my @data = split /:/, $tag;
	  			#print join(':', @data), "\n";
	  			if($data[0] eq "SM") { $sm = $data[1]; }
	  			if($data[0] eq "ID") { $id = $data[1]; }
	  		}
	  		push @{$readgroups->{$sm}->{rg}}, $id;
	  	}
	}
	close HEADER;
}

## save read groups to file
#####################################################
sub saveReadGroup {
	
	my $sample = $_[0];
	my $PREFIX = $_[1];
	my $readgroups = $_[2];
	my $rgfile = $_[3];
	
	open SAMPLE, "> $PREFIX/$rgfile" or die "Can't open $PREFIX/$rgfile ($!)\n";
	
	for my $sm ( keys %$readgroups ) {
		
		if($sample eq "ALL") {
			#print $sample.": ";
			foreach my $rg ( @{$readgroups->{$sm}->{rg}}) {
				#	print $rg.", ";
				print SAMPLE "$rg\n"; 
			}
		}
		elsif($sm eq $sample) {
			#print $sample.": ";
			foreach my $rg ( @{$readgroups->{$sm}->{rg}}) {
				#	print $rg.", ";
				print SAMPLE "$rg\n"; 
			}
		}
	}
	close SAMPLE;
}

1;
