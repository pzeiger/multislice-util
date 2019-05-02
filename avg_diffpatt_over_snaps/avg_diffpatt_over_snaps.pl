#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

# set default output prefixes
my $fileavgdiffout = "avgdiffpatt";


# Initialize arrays and set defaults
my @filein = ();
my @outevery = ();
my ($setin, $setout, $setshift) = (0, 0, 0);
my ($incohout, $cohout) = (1, 1);
#my ($incohout, $cohout, $tdsout) = (1, 1, 1);


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
	if ( $ARGV[$i] =~ /^-noincoh$/ ) {
		$incohout = 0;
		next;
	}
	if ( $ARGV[$i] =~ /^-nocoh$/ ) {
		$cohout = 0;
		next;
	}
#	if ( $ARGV[$i] =~ /^-notds$/ ) {
#		$tdsout = 0;
#		next;
#	}
	
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

if ( !@outevery ) {
	push @outevery, $#filein;
}

my @sumint = ();
my @sumamp = ();
my @Mincoh = ();
my @Sincoh = ();
my @Mcoh = ();
my @Scoh = ();
my $ixmax = 0;
my $iymax = 0;

for my $j ( 0 .. $#filein ) {
	
	my @inp = ();
    my $ixmax2 = 0;
    my $iymax2 = 0;
	
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
				$ixmax2 = $inp[$#inp][1];
			}
			if ( $inp[$#inp][2] >= $iymax ) {
				$iymax2 = $inp[$#inp][2];
			}
		}
	}
    
    if ( not $ixmax ) {
        $ixmax = $ixmax2;
    } elsif ( $ixmax != $ixmax2 ) {
        die;
    }
    if ( not $iymax ) {
        $iymax = $iymax2;
    } elsif ( $iymax != $iymax2 ) {
        die;
    }
    
	foreach my $index ( 0 .. $#inp ) {
		my $ix = $inp[$index][1];
		my $iy = $inp[$index][2];
		my $int;
		my @amp = ();
		
		if ( $setshift == 1 ) {
			$ix = (( $ix - 1 + ($ixmax/2) ) % $ixmax) + 1;
			$iy = (( $iy - 1 + ($iymax/2) ) % $iymax) + 1;
		}
		# incoherent sum
		$int = $inp[$index][3]*$inp[$index][3] + $inp[$index][4]*$inp[$index][4];
		$sumint[$ix][$iy] += $int;
		# compute running stddev according to Knuth's algorithm (https://www.johndcook.com/blog/standard_deviation/)
		if ($j == 0) {
			$Mincoh[$ix][$iy][1] = $int;
			$Sincoh[$ix][$iy][1] = .0;
		} else {
			$Mincoh[$ix][$iy][2] = $Mincoh[$ix][$iy][1] + ($int - $Mincoh[$ix][$iy][1]) / ($j+1);
			$Sincoh[$ix][$iy][2] = $Sincoh[$ix][$iy][1] + ($int - $Mincoh[$ix][$iy][1]) * ($int - $Mincoh[$ix][$iy][2]);
			$Mincoh[$ix][$iy][1] = $Mincoh[$ix][$iy][2];
			$Sincoh[$ix][$iy][1] = $Sincoh[$ix][$iy][2];
		}
		
		# coherent sum
		$amp[1] = $inp[$index][3];
		$amp[2] = $inp[$index][4];
		$sumamp[$ix][$iy][1] += $amp[1];
		$sumamp[$ix][$iy][2] += $amp[2];
		
		# compute running stddev according to Knuth's algorithm (https://www.johndcook.com/blog/standard_deviation/)
		if ($j == 0) {
			$Mcoh[$ix][$iy][1][1] = $amp[1];
			$Scoh[$ix][$iy][1][1] = .0;
			$Mcoh[$ix][$iy][2][1] = $amp[2];
			$Scoh[$ix][$iy][2][1] = .0;
		} else {
			$Mcoh[$ix][$iy][1][2] = $Mcoh[$ix][$iy][1][1] + ($amp[1] - $Mcoh[$ix][$iy][1][1]) / ($j+1);
			$Mcoh[$ix][$iy][2][2] = $Mcoh[$ix][$iy][2][1] + ($amp[2] - $Mcoh[$ix][$iy][2][1]) / ($j+1);

			$Scoh[$ix][$iy][1][2] = $Scoh[$ix][$iy][1][1] + ($amp[1] - $Mcoh[$ix][$iy][1][1]) * ($amp[1] - $Mcoh[$ix][$iy][1][2]);
			$Scoh[$ix][$iy][2][2] = $Scoh[$ix][$iy][2][1] + ($amp[2] - $Mcoh[$ix][$iy][2][1]) * ($amp[2] - $Mcoh[$ix][$iy][2][2]);
			
			$Mcoh[$ix][$iy][1][1] = $Mcoh[$ix][$iy][1][2];
			$Scoh[$ix][$iy][2][1] = $Scoh[$ix][$iy][2][2];
		}
	}
	
	foreach my $n ( @outevery ) {
		if ( $j > 0 and ((( $j + 1 ) % $n) == 0 or $j == $#filein) ) {
			my $outdp;
            
            if ( $incohout or $cohout ) {
                my $FILEOUT = $fileavgdiffout . ($j+1);
                open $outdp, '>:encoding(UTF-8)', $FILEOUT;
#                if ( $j == $#filein ) {
#                    my $FILEOUT = $fileavgdiffout . ($#filein + 1) . "final";
#                    open $outdp, '>:encoding(UTF-8)', $FILEOUT;
#                } else {
#                    my $FILEOUT = $fileavgdiffout . ($j+1);
#                    open $outdp, '>:encoding(UTF-8)', $FILEOUT;
#                }
            }
            
            printf $outdp "# Averaged diffraction pattern. Columns contain:\n";
            printf $outdp "# 1 -> x pixel\n";
            printf $outdp "# 2 -> y pixel\n";
            printf $outdp "# 3 -> incohint\n";
            printf $outdp "# 4 -> uincohint\n";
            printf $outdp "# 5 -> cohint\n";
            printf $outdp "# 6 -> ucohint\n";
            
			my @incohint = ();
			my @uincohint = ();
			my @stddevint = ();
			my @intincohint = ();
            
			my @cohint = ();
			my @ucohint = ();
			my @intcohint = ();
			my @stddevamp = ();
			my @usumamp = ();
			
			my $tdsint;
			my $utdsint;
			my @inttdsint = ();
			
			# calculate incoherent and coherent avg
			for my $ix ( 1 .. $#sumint ) {
				for my $iy ( 1 .. $#{$sumint[$ix]} ) {
					# incoherent avg
					$incohint[$ix][$iy] = $sumint[$ix][$iy] / ($j+1);
					$stddevint[$ix][$iy] = sqrt($Sincoh[$ix][$iy][1]/$j);
					$uincohint[$ix][$iy] = sqrt(1/($j+1)) * $stddevint[$ix][$iy];
					# coherent avg
					$stddevamp[1] = sqrt($Scoh[$ix][$iy][1][1]/$j);
					$stddevamp[2] = sqrt($Scoh[$ix][$iy][2][1]/$j);
					$usumamp[1] = sqrt(1/($j+1)) * $stddevamp[1];
					$usumamp[2] = sqrt(1/($j+1)) * $stddevamp[2];
						
                    
                    if ( $incohout ) {
                        printf $outdp "%4i %4i %1.12e %1.12e ", ($ix, $iy, $incohint[$ix][$iy], $uincohint[$ix][$iy]);
                    } else {
                        printf $outdp "%4i %4i # # ", ($ix, $iy);
                    }
                    
                    if ( $cohout ) {
					    $cohint[$ix][$iy] = ($sumamp[$ix][$iy][1]*$sumamp[$ix][$iy][1] + $sumamp[$ix][$iy][2]*$sumamp[$ix][$iy][2]) / (($j+1)**2);
					    $ucohint[$ix][$iy] = sqrt(( 2*$sumamp[$ix][$iy][1]*$usumamp[1] )**2 + ( 2*$sumamp[$ix][$iy][2]*$usumamp[2] )**2) / (($j+1)**2);
                        printf $outdp "%1.12e %1.12e \n", ($cohint[$ix][$iy], $ucohint[$ix][$iy]);
                    } else {
                        printf $outdp "# # \n";
                    }
                    
#                    if ( $tdsout ) {
#					    $tdsint = $incohint[$ix][$iy] - $cohint[$ix][$iy];
#					    $utdsint = sqrt($uincohint[$ix][$iy]**2 + $ucohint[$ix][$iy]**2);
#                        printf $outdp "%1.12e %1.12e\n", ($tdsint, $utdsint);
#                    } else {
#                        printf $outdp "# #\n";
#                    }
				}
                printf $outdp "\n";
			}
			close ($outdp);
		}
	}
}

