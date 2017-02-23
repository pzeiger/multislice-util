#!/usr/bin/perl


# Calculates the msd of snapshots at given time steps from lammps dump file
# and outputs timestep vs msd 


use strict;
use warnings;

# Number of lines containing environment information per time step in lammps output
my $nlmphead = 8;


my $filein = $ARGV[0];
my $qatom = $ARGV[1];

# turn the center of mass (com) translation correction on (1) or off (all other values)
my $trlcorr = $ARGV[2];
if ( $trlcorr == 1 ) {
	print "COM translation correction turned on\n";
}
else {
	$trlcorr = 0;
	print "COM translation correction turned off\n";
}

# Enable COM correction in PIMD simulations
my $mode = $ARGV[3];
my $centroidfile;
if ( $mode =~ /pimd/ ) {
	print "Operating in PIMD mode\n";
	$centroidfile = $ARGV[4];
}

open my $in, '<', $filein;
my $centroidin;
if ( $mode =~ /pimd/ ) {
	open $centroidin, '<', $centroidfile;
}


# create initial structure array; if given use first snapshot of the
# input file if not extra file ($ARGV[2]) is specified 

my $line = <$in>;
my ($tstep0, $natom0, $coord0, $snap0, $dim0, $ixyz0) = readTimestepData($in,$line,$nlmphead);

my @coord0 = @$coord0;
my @snap0 = @$snap0;
my @dim0 = @$dim0;
my @ixyz0 = @$ixyz0;


my @centroid = ();
my $centroid;

if ( $mode =~ /pimd/ && $trlcorr == 1 ) {
	my $line = <$centroidin>;
	(my $tstep, my $natom, my $coord, $centroid, my $dim, my $ixyz) = readTimestepData($centroidin,$line,$nlmphead);
	
	@centroid = @$centroid;
	
	print "@{$centroid[1]}\n";
	print "$centroid\n";
}

# determine center of mass of the crystal structure
my @com0;
if ( $mode !~ /pimd/ ) {
	@com0 = determineCOM($snap0,$natom0,$ixyz0);
}
else {
	@com0 = determineCOM($centroid,$natom0,$ixyz0);
}

print "CENTRE OF MASS: @com0\n";
print "COORDINATES: @coord0\n";

if ( $mode !~ /pimd/ ) {
	my $fileout = "snapshot" . $tstep0;
	open my $out,'>', $fileout;
	
	printf $out "%f %f %f\n", @dim0;
	printf $out "%i F\n", $natom0;
	foreach ( @snap0 ) {
		printf $out "$qatom @{$_}[$ixyz0[0]] @{$_}[$ixyz0[1]] @{$_}[$ixyz0[2]] \n";
	}
}
else {
	my $fileout = "snapshot" . $tstep0;
	open my $out,'>', $fileout;
	
	printf $out "%f %f %f\n", @dim0;
	printf $out "%i F\n", $natom0;
	foreach ( @centroid ) {
		printf $out "$qatom @{$_}[$ixyz0[0]] @{$_}[$ixyz0[1]] @{$_}[$ixyz0[2]] \n";
	}
}


# Print all other snapshots to multislice input files
while ( my $line = <$in> ) {
	my ($tstep, $natom, $coord, $snap, $dim, $ixyz) = readTimestepData($in,$line,$nlmphead);
	
	my @coord = @$coord;
	my @snap = @$snap;
	my @dim = @$dim;
	
	my @ixyz = @$ixyz;
	
	my @comcorr = ();

	# correct for atoms reentering simulation box on the opposite site
	# due to periodic boundary conditions
	foreach my $i ( 0 .. $#snap ) {
		foreach my $j ( $ixyz[0] .. $ixyz[2] ) {
#			print "$i $j $snap[$i][$j]\n";
#			print "$j @{$snap[$i]}\n";
			my $disp = $snap[$i][$j] - $snap0[$i][$j];
			if ($disp > 0.5) {
				$snap[$i][$j] -= 1;
			}
			if ($disp < - 0.5) {
				$snap[$i][$j] += 1;
			}
		}
	}
	
	
	# In PIMD mode read centroid positions
	if ( $mode =~ /pimd/ ) {
		my $line = <$centroidin>;
		(my $tstep, my $natom, my $coord, $centroid, my $dim, my $ixyz) = readTimestepData($centroidin,$line,$nlmphead);
		@centroid = @$centroid;
		
		print "@{$centroid[1]}\n";
	}
	
	# Calculate necessary center of mass correction
	my @com;
	if ( $mode !~ /pimd/ ) {
		    @com = determineCOM($snap,$natom,$ixyz);
		}
		else {
		    @com = determineCOM($centroid,$natom,$ixyz);
	}

	if ($trlcorr == 1) {
		foreach my $d ( 0 .. $#com ) {
			push @comcorr, ( $com[$d] - $com0[$d] );
		}
	}
	else {
		foreach my $d ( 0 .. $#com ) {
			push @comcorr, 0;
		}
	}
	
	print "\nTIMESTEP: $tstep\n";
	print "CENTRE OF MASS: @com\n";
	print "COMCORR: @comcorr\n";
	
	my $fileout = "snapshot" . $tstep;
	open my $out,'>', $fileout;
	
	printf $out "%f %f %f\n", @dim;
	printf $out "%i F\n", $natom;
	foreach ( @snap ) {
		@{$_}[$ixyz[0]] -= $comcorr[0];
		@{$_}[$ixyz[1]] -= $comcorr[1];
		@{$_}[$ixyz[2]] -= $comcorr[2];
		
		if ( @{$_}[$ixyz[0]] < 0 ) {@{$_}[$ixyz[0]] +=1;}
		if ( @{$_}[$ixyz[1]] < 0 ) {@{$_}[$ixyz[1]] +=1;}
		if ( @{$_}[$ixyz[2]] < 0 ) {@{$_}[$ixyz[2]] +=1;}
		
		if ( @{$_}[$ixyz[0]] > 1 ) {@{$_}[$ixyz[0]] -=1;}
        if ( @{$_}[$ixyz[1]] > 1 ) {@{$_}[$ixyz[1]] -=1;}
        if ( @{$_}[$ixyz[2]] > 1 ) {@{$_}[$ixyz[2]] -=1;}
		
		printf $out "$qatom @{$_}[$ixyz[0]] @{$_}[$ixyz[1]] @{$_}[$ixyz[2]] \n";
	}
}


# subroutine that reads data for 1 time step from a LAMMPS dump file formatted as:
# ------BEGIN FILE------
# ITEM: TIMESTEP
# 0
# ITEM: NUMBER OF ATOMS
# 34680
# 0.0  62.426 xlo xhi
# 0.0  61.271 ylo yhi
# 0.0 300.165 zlo zhi
# ITEM: ATOMS id type xs ys zs
# 1 1 0.0000000 0.0000000 0.0000000
# 2 1 0.0333333 0.0000000 0.0098039
# ------END FILE ------
# 
# takes 3 arguments. Current input line $line, input file identifier $in and number
# of lines $nhead composing the header

sub readTimestepData {
	my $in = $_[0];
	my $line = $_[1];
	my $nlmphead = $_[2];
	
	chomp($line);
	my @head = ($line);
	
	foreach (1 .. $nlmphead) {
		$line = <$in>;
		chomp($line);
		push @head, $line;
	}
	
	# 2nd line of the header contains time step info;
	# 3rd line the number of atoms
	my $tstep = $head[1];
	my $natom = $head[3];
	# arrays for unit cell dimensions, atom data and temporary auxiliary purposes
	my @dim	= ();
	my @snap = ();
	my %data;
	
	my @ixyz = ();
	my $id;
	my @coord;
	
	my @tmp = split ' ', $head[$nlmphead];
	
	foreach my $i ( 0 .. $#tmp ) {
		if ( $tmp[$i] =~ /id/ ) {
			$id = $i;
		}
	}
	
	foreach my $i ( 0 .. $#tmp ) {
		if ( $tmp[$i] =~ /xs/ ) {
			$ixyz[0] = $i - $id;
			$coord[0] = "F";
		}
		if ( $tmp[$i] =~ /^x$/ ) {
			$ixyz[0] = $i - $id;
			$coord[0] = "R";
		}
		if ( $tmp[$i] =~ /ys/ ) {
			$ixyz[1] = $i - $id;
			$coord[1] = "F";
		}
		if ( $tmp[$i] =~ /^y$/ ) {
			$ixyz[1] = $i - $id;
			$coord[1] = "R";
		}
		if ( $tmp[$i] =~ /zs/ ) {
			$ixyz[2] = $i - $id;
			$coord[2] = "F";
		}
		if ( $tmp[$i] =~ /^z$/ ) {
			$ixyz[2] = $i - $id;
			$coord[2] = "R";
		}
	}
	
	foreach my $i ( ($nlmphead-3) .. ($nlmphead-1) ) {
		@tmp = split ' ', $head[$i];
		push @dim, ( $tmp[1] - $tmp[0] );
	}
	
	foreach ( 0 .. $natom - 1 ) {
		$line = <$in>;
		chomp($line);
		@tmp = split ' ', $line;
		
		$data{$tmp[0]} = $line;
	}
	
	foreach my $id (1 .. $natom) {
		@tmp = split ' ', $data{$id};
		push @snap, [ @tmp ];
	}
	
	# \@ creates a reference to the array -> see:
	# http://perlmeme.org/faqs/perl_thinking/returning.html
	return($tstep, $natom, \@coord, \@snap, \@dim, \@ixyz);
}


# determine center of mass of a given snapshot

sub determineCOM {
	# first argument contains reference name to array
	my $snap = $_[0];
	my $natom = $_[1];
	my $ixyz = $_[2];
	
	my @snap = @$snap;
	my @ixyz = @$ixyz;
	
	my $comx;
	my $comy;
	my $comz;
	
	foreach my $id ( 0 .. $#snap) {
		$comx += $snap[$id][$ixyz[0]];
		$comy += $snap[$id][$ixyz[1]];
		$comz += $snap[$id][$ixyz[2]];
	}
	$comx /= $natom;
	$comy /= $natom;
	$comz /= $natom;
	my @com = ($comx, $comy, $comz);
	
	return (@com);
}


