package Utils;

###################################################################
# Utils
#
# Package with basic commons routines
# 
#  Author: Giuseppe Narzisi 
#    Date: December 11, 2013
#
###################################################################

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(runCmd normalizeIndel uniq binarySearch inTarget bychrpos elapsedTime);
@EXPORT_OK = qw($findVariants $findDenovos $findSomatic $exportTool $bamtools);

use strict;
use warnings;
use POSIX;
use FindBin qw($Bin);
use lib $Bin; # add $Bin directory to @INC

# programs via absolute path
our $findVariants = "$Bin/FindVariants.pl";
our $findDenovos  = "$Bin/FindDenovos.pl";
our $findSomatic  = "$Bin/FindSomatic.pl";
our $exportTool   = "$Bin/ExportVariants.pl";
our $bamtools     = "$Bin/bamtools-2.3.0/bin/bamtools";

# Run system command 
#####################################################
sub runCmd 
{ 
	my $info = $_[0];
	my $cmd  = $_[1];
	#my $VERBOSE = $_[2];

	#print STDERR "$info ($cmd)...\n" if $VERBOSE;

	my $rc = system($cmd); 
	#die "failed ($rc)\n" if $rc; 
	if ($rc) { 
		print STDERR "Command failure: $info ($cmd)...\n";

		if ($? == -1) {
			print "failed to execute: $!\n";
		}
		elsif ($? & 127) {
			printf "child died with signal %d, %s coredump\n",
			($? & 127),  ($? & 128) ? 'with' : 'without';
		}
		else {
			printf "child exited with value %d\n", $? >> 8;
		}	
	}
}

## Normalize indel location to leftmost position
#####################################################
sub normalizeIndel {
	my $denovo  = ${$_[0]} ;
	my $delta = $_[1];
	my $direction = $_[2];
	my $genome = $_[3];	
	
	my $pos = $denovo->{pos};
	my $type = $denovo->{type};
	my $chr = $denovo->{chr};

	## extract sequence around indel from reference file
	#$pos = $pos+1; # convert from 0-based to 1-based coordinate system
	my $l = $pos-$delta;
 	my $r = $pos+$delta;
	#print "samtools faidx $REF $chr:$l-$r | awk '\$0 ~ /^>/' \n";
	#my $reference = readpipe( "$samtools faidx $REF $chr:$l-$r | awk '\$0 !~ /^>/'" ); # 0-based coordinate system
	my $reference = substr($genome->{$chr}->{seq}, $l-1, $r-$l+1);
	#print "samtools faidx $REF $chr:$l-$r\n";

 	my $p = $delta;
	my $indel_len = $denovo->{len};
	my $new_pos = $pos;
	if( $type eq "del") { # deletion
		my $left_seq  = substr($reference, 0, $p);
		my $right_seq = substr($reference, $p+$indel_len);
		my $true_indel_seq = $left_seq . $right_seq;
		$new_pos = $l+$p;
		$denovo->{ref} = substr($reference, $p, $indel_len);

		## reposition the candidate variant to a canonical (leftmost/rightmost) position
		my $left_hyplotype;
		my $right_hyplotype;
		my $new_indel_seq;
		for (my $i = $p-1; $i >= 0; $i--) {
			$left_hyplotype  = substr($reference, 0, $i);
		  	$right_hyplotype = substr($reference, $i+$indel_len);
		  	$new_indel_seq = $left_hyplotype . $right_hyplotype;
		  	#print "$true_indel_seq";
		    #print "$new_indel_seq";
		   	if ($true_indel_seq eq $new_indel_seq) {
		   		#print "Repositioning!\n";
		  		$new_pos = $l+$i; # repositioning!
		  		$denovo->{ref} = substr($reference, $i, $indel_len);
		  	}
		}
		$denovo->{seq} = ( '-' x $indel_len );
	}
	# Simulator behaviour: when deleting
	# a string it starts deleting at position i, while, when inserting a
	# string, it appends it starting at position i+1.
	elsif( $type eq "ins") {
	  	my $left_seq  = substr($reference, 0, $p);
		my $right_seq = substr($reference, $p);
		my $true_indel_seq = $left_seq . $denovo->{seq} . $right_seq;
		$new_pos = $l+$p;

		## reposition the candidate variant to a canonical (leftmost/rightmost) position
		my $left_hyplotype;
		my $right_hyplotype;
		my $new_indel_seq;
		my $ins;
		for (my $i = $p-1; $i >= 0; $i--) {
			$left_hyplotype  = substr($reference, 0, $i);
		  	$right_hyplotype = substr($reference, $i);
		  	$ins = substr($true_indel_seq, $i, $indel_len);
		  	$new_indel_seq = $left_hyplotype . $ins . $right_hyplotype;
		  	#print "$true_indel_seq";
		    #print "$new_indel_seq";
		  	if ($true_indel_seq eq $new_indel_seq) {
		  		#print "Repositioning!\n";
		  		$new_pos = $l+$i; # repositioning!
		  		$denovo->{seq} = $ins;
		  	}
		}
		$denovo->{ref} = ( '-' x $indel_len );
	}
	$denovo->{pos} = $new_pos;
}

## return unique elements of input list
#####################################################
sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}

# binary-search
##########################################
sub binarySearch {

	my $x = $_[0];	
	my $A = $_[1];
	my $l = $_[2];
	my $r = $_[3];
	
	while($l<=$r) {
		my $m = floor( ($l+$r)/2 );
		my $f = @$A[$m];
		if($f->{end} < $x) { $l = $m+1; }
		elsif($f->{start} > $x)  { $r = $m-1; }
		else { 
			return $m; 
		}
	}
	return -1;
}

## return true if mutation is in target (start or end
## positions intersect the target), false otherwise
#####################################################
sub inTarget {
	my $mutation = $_[0];
	my $exons = $_[1];
	
	my $chr = $mutation->{chr};
	my $pos = $mutation->{pos};
	my $len = $mutation->{len};
	
	my $result = "false";
	
	my @array = @{$exons->{$chr}};
	
	my $flag1 = binarySearch($pos, \@array, 0, (scalar @array)-1);
	my $flag2 = binarySearch($pos+$len, \@array, 0, (scalar @array)-1);
	
	if( ($flag1 != -1) || ($flag2 != -1) ) { $result = "true"; }
	
	#foreach my $exon (@{$exons->{$chr}}) {
		
	#	my $s = $exon->{start};
	#	my $e = $exon->{end};
		
	#	if( ($s <= $pos) && ($pos <= $e) ) { 
	#		$result = "true"; 
	#		last;
	#	}
	#}
	
	return $result;
}

# auxilary routine to sort variants (Lexicographically)
# by chromosome and then by position
##########################################

sub bychrpos {
	my ($a, $b, $h) = @_;
	my $result;
	if ( ($h->{$a}->{chr} =~ m/^-?\d+\z/) && ($h->{$b}->{chr} =~ m/^-?\d+\z/) ) {
		$result = ( ($h->{$a}->{chr} <=> $h->{$b}->{chr}) || ($h->{$a}->{pos} <=> $h->{$b}->{pos}) );
	}
	else {
		$result = ( ($h->{$a}->{chr} cmp $h->{$b}->{chr}) || ($h->{$a}->{pos} <=> $h->{$b}->{pos}) );
	}
	return $result;
}

# print total time elapsed 
##########################################
sub elapsedTime {
	my $time_taken = $_[0];
	my $tool = $_[1];

	my $hours   = $time_taken / 3600;
	my $days = $hours / 24; 
	my $seconds = $time_taken % 3600;
	my $minutes = $seconds / 60;
	$seconds    = $seconds % 60;

	$days = floor($days);
	$hours = floor($hours);
	$minutes = floor($minutes);
	$seconds = floor($seconds);

	print STDERR "$tool elapsed time: $days day(s) $hours hour(s) $minutes minute(s) $seconds second(s)\n";
}

# log 10
##########################################
#sub log10 {
#        my $n = shift;
#        return log($n)/log(10);
#}

1;