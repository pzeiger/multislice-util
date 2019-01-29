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


my @avgint = ();
#my @avgsqint = ();
my @M = ();
my @S = ();


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
		my $int;
		
		if ( $setshift == 1 ) {
			$ix = (( $ix - 1 + ($ixmax/2) ) % $ixmax) + 1;
			$iy = (( $iy - 1 + ($iymax/2) ) % $iymax) + 1;
		}
		$int = $inp[$index][3]*$inp[$index][3] + $inp[$index][4]*$inp[$index][4];
		$avgint[$ix][$iy] += $int;
		
#		$avgsqint[$ix][$iy] += ($inp[$index][3]*$inp[$index][3] + $inp[$index][4]*$inp[$index][4])**2;
#		$avgsqint[$ix][$iy] += ($inp[$index][3]*$inp[$index][3] + $inp[$index][4]*$inp[$index][4])**2;
		
		# compute running stddev according to Knuth's algorithm (https://www.johndcook.com/blog/standard_deviation/)
		if ($j == 0) {
			$M[$ix][$iy][1] = $int;
			$S[$ix][$iy][1] = .0;
		} else {
			$M[$ix][$iy][2] = $M[$ix][$iy][1] + ($int - $M[$ix][$iy][1]) / ($j+1);
			$S[$ix][$iy][2] = $S[$ix][$iy][1] + ($int - $M[$ix][$iy][1]) * ($int - $M[$ix][$iy][2]);
			$M[$ix][$iy][1] = $M[$ix][$iy][2];
			$S[$ix][$iy][1] = $S[$ix][$iy][2];
		}
	}
	
	foreach my $n ( @outevery ) {
		if ( (( $j + 1 ) % $n) == 0 or $j == $#filein ) {
			my $stddev;
			my $uncert;
			my $out;
#			for my $ix  ( 1 .. $#avgint ) {
#				for my $iy ( 1 .. $#{$avgint[$ix]} ) {
#					$stddev[$ix][$iy] += sqrt(1/$j * abs(($avgsqint[$ix][$iy] - $avgint[$ix][$iy]**2/($j+1))));
#					$uncavgint[$ix][$iy] += sqrt(1/($j+1)) * $stddev[$ix][$iy];
#				}
#			}
			if ( $j == $#filein ) {
				my $FILEOUT = "incohavg_diffpatt" . ($#filein + 1) . "final";
				open $out, '>:encoding(UTF-8)', $FILEOUT;
			} else {
				my $FILEOUT = "incohavg_diffpatt" . ($j+1);
				open $out, '>:encoding(UTF-8)', $FILEOUT;
			}
			
			for my $ix ( 1 .. $#avgint ) {
				for my $iy ( 1 .. $#{$avgint[$ix]} ) {
					$stddev = sqrt($S[$ix][$iy][1]/$j);
					$uncert = sqrt(1/($j+1)) * $stddev;
					printf $out "%4i %4i %1.12e %1.12e %1.12e\n", ($ix, $iy, $avgint[$ix][$iy]/($j+1), $stddev, $uncert);
				}
				printf $out "\n";
			}
			close ($out);
		}
	}
}

