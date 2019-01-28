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
my @int = ();


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
			$avgint[$ix][$iy][1] += $inp[$index][3]
			$avgint[$ix][$iy][2] += $inp[$index][4]
#			$int[$j][$ix][$iy] = $inp[$index][3]*$inp[$index][3] + $inp[$index][4]*$inp[$index][4];
	}
	
	foreach my $n ( @outevery ) {
		if ( (( $j + 1 ) % $n) == 0 ) {
#			my @stddev = ();
#			for my $i ( 0 .. $j ) {
#				for my $ix  ( 1 .. $#avgint ) {
#					for my $iy ( 1 .. $#{$avgint[$ix]} ) {
#						$stddev[$ix][$iy] += ($int[$i][$ix][$iy] - ($avgint[$ix][$iy]/($j+1)))**2;
#					}
#				}
#			}
			my $FILEOUT = "cohavg_diffpatt" . ($j+1);
			open my $out, '>:encoding(UTF-8)', $FILEOUT;
			
			for my $ix ( 1 .. $#avgint ) {
				for my $iy ( 1 .. $#{$avgint[$ix]} ) {
					printf $out "%4i %4i %1.12e %1.12e %1.12e\n", ($ix, $iy, ($avgint[$ix][$iy][1]*$avgint[$ix][$iy][1] + $avgint[$ix][$iy][2]*$avgint[$ix][$iy][2])/($j+1));
				}
				printf $out "\n";
			}
			close ($out);
		}
	}
}


#my @stddev = ();
#
#for my $j ( 0 .. $#filein ) {
#	for my $ix  ( 1 .. $#avgint ) {
#		for my $iy ( 1 .. $#{$avgint[$ix]} ) {
#			$stddev[$ix][$iy] += ($int[$j][$ix][$iy] - ($avgint[$ix][$iy]/($#filein+1)))**2;
#		}
#	}
#}

my $FILEOUT = "avg_diffpatt" . ($#filein + 1) . "final";

open my $out, '>:encoding(UTF-8)', $FILEOUT;

for my $ix ( 1 .. $#avgint ) {
	for my $iy ( 1 .. $#{$avgint[$ix]} ) {
		printf $out "%4i %4i %1.12e %1.12e %1.12e\n", ($ix, $iy, ($avgint[$ix][$iy][1]*$avgint[$ix][$iy][1] + $avgint[$ix][$iy][2]*$avgint[$ix][$iy][2])/($j+1));#, (sqrt($stddev[$ix][$iy]/($#filein))), (sqrt($stddev[$ix][$iy]/(($#filein)*($#filein+1)))));
	}
	printf $out "\n";
}
close ($out);


