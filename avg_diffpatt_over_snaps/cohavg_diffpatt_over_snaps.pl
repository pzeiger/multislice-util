#!/usr/bin/perl


use strict;
use warnings;


# Initialize arrays and set defaults
my @filein = ();
my @outevery = ();
my ($setin, $setout, $setshift) = (0, 0, 0);


# Loop over arguments and place them into arrays for further use
foreach my $i ( 0 .. $#ARGV ) {

	# Recognize which settings are to be made
	if ( $ARGV[$i] =~ /^-in$/ ) {
		($setin, $setout) = (1, 0);
		next;
	}
	if ( $ARGV[$i] =~ /^-outevery$/ ) {
		($setin, $setout) = (0, 1,);
		next;
	}
	if ( $ARGV[$i] =~ /^-shift$/ ) {
		$setshift = 1;
		next;
	}
	
	# set all the calculation parameters
	if ( $setin == 1 ) {
		push @filein, $ARGV[$i];
		next;
	}
	if ( $setout == 1 ) {
		push @outevery, $ARGV[$i];
		next;
	}
	
}


my @sumamp = ();
my @sumsqamp = ();
my $int;
my $uncint;


for my $j ( 0 .. $#filein ) {
	
	my @inp = ();
	my $ixmax = 0;
	my $iymax = 0;
	
	open my $in, '<', $filein[$j]
		or do{
			warn "Could not open file '$filein[$j]' $!";
			next;
		};
	while ( <$in> ) {
		if ( $_ =~ /^\s$/ ) {
			next;
		}
		else {
			push @inp, [ split /\s+/ ];
			if ( $inp[$#inp][1] >= $ixmax ) {
				$ixmax = $inp[$#inp][1];
			}
			if ( $inp[$#inp][2] >= $iymax ) {
				$iymax = $inp[$#inp][2];
			}	
		}
	}
	foreach my $index ( 0 .. $#inp ) {
		my $ix = $inp[$index][1];
		my $iy = $inp[$index][2];
		
		if ( $setshift == 1 ) {
			$ix = (( $ix - 1 + ($ixmax/2) ) % $ixmax) + 1;
			$iy = (( $iy - 1 + ($iymax/2) ) % $iymax) + 1;
		}
			$sumamp[$ix][$iy][1] += $inp[$index][3];
			$sumamp[$ix][$iy][2] += $inp[$index][4];
			$sumsqamp[$ix][$iy][1] += $inp[$index][3]*$inp[$index][3];
			$sumsqamp[$ix][$iy][2] += $inp[$index][4]*$inp[$index][4];
	}
	
	foreach my $n ( @outevery ) {
		if ( (( $j + 1 ) % $n) == 0 or $j == $#filein ) {
			my @uncert = ();
			my $out;
#			for my $i ( 0 .. $j ) {
			for my $ix  ( 1 .. $#sumamp ) {
				for my $iy ( 1 .. $#{$sumamp[$ix]} ) {
					$uncert[$ix][$iy][1] += sqrt(1/($j**2-$j)*($sumsqamp[$ix][$iy][1] - ($sumamp[$ix][$iy][1])**2/($j+1)));
					$uncert[$ix][$iy][2] += sqrt(1/($j**2-$j)*($sumsqamp[$ix][$iy][2] - ($sumamp[$ix][$iy][2])**2/($j+1)));
				}
			}
#			}
			if ( $j == $#filein ) {
				my $FILEOUT = "cohavg_diffpatt" . ($#filein + 1) . "final";
				open $out, '>:encoding(UTF-8)', $FILEOUT;
			} else {
				my $FILEOUT = "cohavg_diffpatt" . ($j+1);
				open $out, '>:encoding(UTF-8)', $FILEOUT;
			}
			
			for my $ix ( 1 .. $#sumamp ) {
				for my $iy ( 1 .. $#{$sumamp[$ix]} ) {
					$int = ($sumamp[$ix][$iy][1]*$sumamp[$ix][$iy][1] + $sumamp[$ix][$iy][2]*$sumamp[$ix][$iy][2])/($j+1);
					$uncint = sqrt(( 2*$sumamp[$ix][$iy][1]*$uncert[$ix][$iy][1] )**2 + ( 2*$sumamp[$ix][$iy][2]*$uncert[$ix][$iy][2] )**2) / ($j+1);
					printf $out "%4i %4i %1.12e %1.12e\n", ($ix, $iy, $int, $uncint);
				}
				printf $out "\n";
			}
			close ($out);
		}
	}
}

