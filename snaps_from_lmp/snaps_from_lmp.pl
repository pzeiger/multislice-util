#!/usr/bin/perl


# Calculates the msd of snapshots at given time steps from lammps dump file
# and outputs timestep vs msd 


use strict;
use warnings;

# Load command inputs
my $filein = $ARGV[0];
my $qatom = $ARGV[1];

open my $in, '<', $filein;

# Read all lammps atom data output for a time step & convert to multislice snapshot 
while ( my $line = <$in> ) {
	my ($tstep, $natom, $snap, $dim) = readTimestepData($in,$line);
	
	my @snap = @$snap;
	my @dim = @$dim;
	
	my $fileout = "snapshot" . $tstep;
	open my $out,'>', $fileout;
	
	# correct positions for atoms which reenter the simulation box on the opposite site
	# due to periodic boundary conditions such that all atoms are within box
	foreach my $i ( 0 .. $#snap ) {
		foreach my $j ( 2 .. $#{$snap[$i]} ) {
			if ( $snap[$i][$j] < 0 ) {
				$snap[$i][$j] += 1;
			}
			if ( $snap[$i][$j] > 1 ) {
				$snap[$i][$j] -= 1;
			}
		}
	}
	printf $out "%f %f %f\n", @dim;
	printf $out "%i F\n", $natom;
	foreach ( @snap ) {
		printf $out "$qatom @{$_}[2] @{$_}[3] @{$_}[4] \n";
	}
	die;
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
	
	chomp($line);
	my @head = ($line);
	
	foreach (1 .. 8) {
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
	my @tmp = ();
	my %data;
	
	foreach my $i ( 5 .. 7 ) {
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
	return($tstep, $natom, \@snap, \@dim);
}


