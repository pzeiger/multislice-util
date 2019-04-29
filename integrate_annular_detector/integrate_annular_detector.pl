#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;

# Natural constants
use constant {
    PI => 4 * atan2(1, 1),                  # pi
    ELEMENTARY_CHARGE => 1.6021766208e-19,  # elementary charge in C
    MASS_ELECTRON => 9.10938356e-31,        # mass electron in kg
    PLANCK_CONSTANT => 6.62607015e-34,      # Planck constant in Js
    SPEED_OF_LIGHT => 299792458,            # vacuum speed of light in m/s
};


# set default output prefixes
my $intintout = "intdetint";


# Initialize arrays and set defaults
my @filein = ();
my $fileintout = "diffpatt_integrated";
my ($setin, $setout, $setshift, $radprofout) = (0, 0, 0, 0);
my ($incohout, $cohout, $tdsout) = (1, 1, 1);
my ($innerdetang, $outerdetang) = (0.0, 0.0);
my ($detcangx, $detcangy) = 0.0;
my ($ly, $lx) = (1.0, 1.0); # in angstroms
my ($nx, $ny) = (0, 0);
my $accv = 100e3; # in V
my @tmp = ();

# Loop over arguments and place them into arrays for further use
foreach my $i ( 0 .. $#ARGV ) {
    
	if ( $ARGV[$i] =~ /^-/ ) {
		$setin = 0;
	}
    
	# Recognize which settings are to be made
	if ( $ARGV[$i] =~ /^-in$/ ) {
		$setin = 1;
		next;
	}
	# Recognize which settings are to be made
	if ( $ARGV[$i] =~ /^-out$/ ) {
		$fileintout = $ARGV[$i+1];
		next;
	}
	if ( $ARGV[$i] =~ /^-shift$/ ) {
		$setshift = 1;
		next;
	}
	if ( $ARGV[$i] =~ /^-oang$/ ) {
		$outerdetang = $ARGV[$i+1];
		next;
	}
	if ( $ARGV[$i] =~ /^-iang$/ ) {
		$innerdetang = $ARGV[$i+1];
		next;
	}
	if ( $ARGV[$i] =~ /^-detcangx$/ ) {
		$detcangx = $ARGV[$i+1];
		next;
	}
	if ( $ARGV[$i] =~ /^-detcangy$/ ) {
		$detcangy = $ARGV[$i+1];
		next;
	}
    # length simcell in x-dir in angstroms
	if ( $ARGV[$i] =~ /^-lx$/ ) {
		$lx = $ARGV[$i+1];
		next;
	}
    # length simcell in y-dir in angstroms
	if ( $ARGV[$i] =~ /^-ly$/ ) {
		$ly = $ARGV[$i+1];
		next;
	}
    # numerical grid used for multislice
	if ( $ARGV[$i] =~ /^-grid$/ ) {
        print "$ARGV[$i+1]\n";
        @tmp = split('x', $ARGV[$i+1]);
		$nx = $tmp[0];
		$ny = $tmp[1];
		next;
	}
    # get acceleration voltage in V
	if ( $ARGV[$i] =~ /^-accv$/ ) {
		$accv = $ARGV[$i+1];
		next;
	}
	
	# set all the calculation parameters
	if ( $setin == 1 ) {
		push @filein, $ARGV[$i];
		next;
	}
}

print "$lx, $ly, $accv, $nx, $ny, $detcangx, $detcangy, $innerdetang, $outerdetang\n";

# Calculate wave length and scattering angles
my $csq = SPEED_OF_LIGHT * SPEED_OF_LIGHT;
my $r = MASS_ELECTRON * $csq;                                       # rest energy
my $t = $accv * ELEMENTARY_CHARGE;                            # kinetic energy, accv in keV
my $v = SPEED_OF_LIGHT * sqrt(1-(($r*$r)/(($r+$t)*($r+$t))));       # relativistic velocity
print $v / SPEED_OF_LIGHT;
my $mr = MASS_ELECTRON / sqrt(1-(($v*$v)/($csq)));                 # relativistic mass
my $lamb = 1e10 * PLANCK_CONSTANT / ($mr*$v);                       # wave length in angstroms
print "\n$lamb\n";
my $k = PI / $lamb;                                             # length of k-vector
my $dkx = PI / ($lx);                                   # thetax-spacing in mrad
my $dky = PI / ($ly);                                   # thetay-spacing in mrad

#my $detckx = $detcangx*$k;  # $detcangx in mrad
#my $detcky = $detcangy*$k;  # $detcangy in mrad

#print "$detckx, $detcky\n";

my $dpcnx = int($nx/2)+1;
my $dpcny = int($ny/2)+1;
print "$dpcnx, $dpcny\n";


open my $intout, '>', $fileintout
    or do{
        warn "Could not open file ${fileintout} $!";
        next;
    };

my $sumincoh = 0.0;
my $usumincoh = 0.0;
my $sumcoh = 0.0;
my $usumcoh = 0.0;
my $sumtdsint = 0.0;
my $usumtdsint = 0.0;
my @header = ();

my $count = 0;

for my $j ( 0 .. $#filein ) {
	
	my @inp = ();
    my $ixmax2 = 0;
    my $iymax2 = 0;
    
    print "opening $filein[$j] \n";
	open my $in, '<', $filein[$j]
		or do{
			warn "Could not open file $filein[$j] $!";
			next;
		};
    
    my $roifileout = $filein[$j] . '_detc' . $detcangx . '-' . $detcangy . '_cang' . $innerdetang . '-'. $outerdetang . '_roi';
    
	open my $roiout, '>', $roifileout
        or do{
			warn "Could not open file ${roifileout} $!";
			next;
		};
	while ( <$in> ) {
		if ( $_ =~ /^\s$/ ) {
			next;
		} elsif ( $_ =~ /^\#/) {
            push @header, $_;
            print "$_";
            next;
		} else {
			@inp = grep { /\S/ } split(/\s+/, $_);
#			if ( $inp[$#inp][1] >= $ixmax ) {
#				$ixmax2 = $inp[$#inp][1];
#			}
#			if ( $inp[$#inp][2] >= $iymax ) {
#				$iymax2 = $inp[$#inp][2];
#			}
            
            my $thetax = ($inp[0]-$dpcnx)*$dkx/$k*1e3;  # in mrad
            my $thetax2 = asin(($inp[0]-$dpcnx)*$dkx/$k)*1e3;  # in mrad
#            print $thetax, $detcangx, $innerdetang, "\n";
            
            if ( abs($thetax-$detcangx) >= $innerdetang and abs($thetax-$detcangx) <= $outerdetang ) {
                my $thetay = ($inp[1]-$dpcny)*$dky/$k*1e3;  # in mrad
                if ( abs($thetay-$detcangy) >= $innerdetang and abs($thetay-$detcangy) <= $outerdetang ) {
#                    print "$thetax $detcangx $innerdetang \n";
#                    print "$thetay $detcangy $innerdetang \n";
#                    print "quality of small angle approx ", abs($thetax-$thetax2), "\n";
                    my $thetadet = sqrt(($thetax-$detcangx)**2 + ($thetay-$detcangy)**2);
#                    print "$thetadet\n";
                    if ( $thetadet >= $innerdetang and $thetadet <= $outerdetang ) {
#                        print "$count ok \n";
                        $count += 1;
                        $sumincoh += $inp[2];
                        $usumincoh += $inp[3]*$inp[3];
                        $sumcoh += $inp[4];
                        $usumcoh += $inp[5]*$inp[5];
                        $sumtdsint += $inp[6];
                        $usumtdsint += $inp[7]*$inp[7];
                    }
                    printf $roiout "%4i %4i %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e\n", $inp[0], $inp[1], $inp[2], $inp[3], $inp[4], $inp[5], $inp[6], $inp[7];
#                    print "$inp[0], $inp[1], $inp[2], $inp[3], $inp[4], $inp[5], $inp[6], $inp[7]";
                }
            }
		}
	}
    close($roiout);
    
    my $beamposx;
    my $beamposy;
    my @tmp = split('/', $filein[$j]);
    for my $i (0 .. $#tmp) {
        if ( $tmp[$i] =~ /^beampos_/ ) {
            my @tmp2 = split(/_x|_y/, $tmp[$i]);
            $beamposx = $tmp2[1];
            $beamposy = $tmp2[2];
        }
    }
    print "$beamposx, $beamposy\n";
    printf $intout "%2i %2i %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e\n", $beamposx, $beamposy, $sumincoh, sqrt($usumincoh), $sumcoh, sqrt($usumcoh), $sumtdsint, sqrt($usumtdsint); 
}


#    if ( not $ixmax ) {
#        $ixmax = $ixmax2;
#    } elsif ( $ixmax != $ixmax2 ) {
#        die;
#    }
#    if ( not $iymax ) {
#        $iymax = $iymax2;
#    } elsif ( $iymax != $iymax2 ) {
#        die;
#    }
    
	
#	foreach my $n ( @outevery ) {
#		if ( $j > 0 and ((( $j + 1 ) % $n) == 0 or $j == $#filein) ) {
#            
#            if ( $incohout or $cohout or $tdsout) {
#                if ( $j == $#filein ) {
#                    my $FILEOUT = $fileavgdiffout . ($#filein + 1) . "final";
#                    open $outdp, '>:encoding(UTF-8)', $FILEOUT;
#                } else {
#                    my $FILEOUT = $fileavgdiffout . ($j+1);
#                    open $outdp, '>:encoding(UTF-8)', $FILEOUT;
#                }
#            }
#            
#			my @incohint = ();
#			my @uincohint = ();
#			my @stddevint = ();
#			my @intincohint = ();
#            
#			my @cohint = ();
#			my @ucohint = ();
#			my @intcohint = ();
#			my @stddevamp = ();
#			my @usumamp = ();
#			
#			my $tdsint;
#			my $utdsint;
#			my @inttdsint = ();
#			
#			my $out;
#			# calculate incoherent and coherent avg
#			for my $ix ( 1 .. $#sumint ) {
#				for my $iy ( 1 .. $#{$sumint[$ix]} ) {
#					# incoherent avg
#					$incohint[$ix][$iy] = $sumint[$ix][$iy] / ($j+1);
#					$stddevint[$ix][$iy] = sqrt($Sincoh[$ix][$iy][1]/$j);
#					$uincohint[$ix][$iy] = sqrt(1/($j+1)) * $stddevint[$ix][$iy];
#					# coherent avg
#					$stddevamp[1] = sqrt($Scoh[$ix][$iy][1][1]/$j);
#					$stddevamp[2] = sqrt($Scoh[$ix][$iy][2][1]/$j);
#					$usumamp[1] = sqrt(1/($j+1)) * $stddevamp[1];
#					$usumamp[2] = sqrt(1/($j+1)) * $stddevamp[2];
#						
#					$cohint[$ix][$iy] = ($sumamp[$ix][$iy][1]*$sumamp[$ix][$iy][1] + $sumamp[$ix][$iy][2]*$sumamp[$ix][$iy][2]) / (($j+1)**2);
#					$ucohint[$ix][$iy] = sqrt(( 2*$sumamp[$ix][$iy][1]*$usumamp[1] )**2 + ( 2*$sumamp[$ix][$iy][2]*$usumamp[2] )**2) / (($j+1)**2);
#					$tdsint = $incohint[$ix][$iy] - $cohint[$ix][$iy];
#					$utdsint = sqrt($uincohint[$ix][$iy]**2 + $ucohint[$ix][$iy]**2);
#                    
#                    if ( $incohout or $cohout or $tdsout) {
#                        if ( $incohout ) {
#                            printf $outdp "%4i %4i %1.12e %1.12e ", ($ix, $iy, $incohint[$ix][$iy], $uincohint[$ix][$iy]);
#                        } else {
#                            printf $outdp "%4i %4i # # ", ($ix, $iy);
#                        }
#                        
#                        if ( $cohout ) {
#                            printf $outdp "%1.12e %1.12e ", ($cohint[$ix][$iy], $ucohint[$ix][$iy]);
#                        } else {
#                            printf $outdp "# # ";
#                        }
#                        
#                        if ( $tdsout ) {
#                            printf $outdp "%1.12e %1.12e\n", ($tdsint, $utdsint);
#                        } else {
#                            printf $outdp "# #\n";
#                        }
#                    }
#                    
#                    if ( $detout ) {
#                        my $theta = (($ix*$dthetax)-$detcangx)**2 + (($iy*$dthetay)-$detcangy)**2
#                        if ( $theta >= $innerdetang and $theta <= $outerdetang ) {
#                             
#                        }
#                    }
#				}
#			}
#		close ($outdp);
#		}
#	}
#
#			# print incohavg
#			if ( $incohout ) {
#				if ( $j == $#filein ) {
#					my $FILEOUT = $fileoutincoh . ($#filein + 1) . "final";
#					open $out, '>:encoding(UTF-8)', $FILEOUT;
#				} else {
#					my $FILEOUT = $fileoutincoh . ($j+1);
#					open $out, '>:encoding(UTF-8)', $FILEOUT;
#				}
#				for my $ix ( 1 .. $#sumint ) {
#					for my $iy ( 1 .. $#{$sumint[$ix]} ) {
#						printf $out "%4i %4i %1.12e %1.12e %1.12e\n", ($ix, $iy, $incohint[$ix][$iy], $stddevint[$ix][$iy], $uincohint[$ix][$iy]);
#					}
#					printf $out "\n";
#				}
#				close ($out);
#			}
#			
#			# print cohavg
#			if ( $cohout ) {
#				if ( $j == $#filein ) {
#					my $FILEOUT = $fileoutcoh . ($#filein + 1) . "final";
#					open $out, '>:encoding(UTF-8)', $FILEOUT;
#				} else {
#					my $FILEOUT = $fileoutcoh . ($j+1);
#					open $out, '>:encoding(UTF-8)', $FILEOUT;
#				}
#				
#				for my $ix ( 1 .. $#sumamp ) {
#					for my $iy ( 1 .. $#{$sumamp[$ix]} ) {
#						printf $out "%4i %4i %1.12e %1.12e\n", ($ix, $iy, $cohint[$ix][$iy], $ucohint[$ix][$iy]);
#					}
#					printf $out "\n";
#				}
#				close ($out);
#			}
#			
#			# print tdsavg
#			if ( $tdsout ) {
#				if ( $j == $#filein ) {
#					my $FILEOUT = $fileouttds . ($#filein + 1) . "final";
#					open $out, '>:encoding(UTF-8)', $FILEOUT;
#				} else {
#					my $FILEOUT = $fileouttds . ($j+1);
#					open $out, '>:encoding(UTF-8)', $FILEOUT;
#				}
#				
#				for my $ix ( 1 .. $#incohint ) {
#					for my $iy ( 1 .. $#{$incohint[$ix]} ) {
#						$tdsint = $incohint[$ix][$iy] - $cohint[$ix][$iy];
#						$utdsint = sqrt($uincohint[$ix][$iy]**2 + $ucohint[$ix][$iy]**2);
#						
#						printf $out "%4i %4i %1.12e %1.12e\n", ($ix, $iy, $tdsint, $utdsint);
#					}
#					printf $out "\n";
#				}
#				close ($out);
#			}
#
#
#
#
#	foreach my $index ( 0 .. $#inp ) {
#		my $ix = $inp[$index][1];
#		my $iy = $inp[$index][2];
#		my $int;
#		my @amp = ();
#		
#		if ( $setshift == 1 ) {
#			$ix = (( $ix - 1 + ($ixmax/2) ) % $ixmax) + 1;
#			$iy = (( $iy - 1 + ($iymax/2) ) % $iymax) + 1;
#		}
#		# incoherent sum
#		$int = $inp[$index][3]*$inp[$index][3] + $inp[$index][4]*$inp[$index][4];
#		$sumint[$ix][$iy] += $int;
#		# compute running stddev according to Knuth's algorithm (https://www.johndcook.com/blog/standard_deviation/)
#		if ($j == 0) {
#			$Mincoh[$ix][$iy][1] = $int;
#			$Sincoh[$ix][$iy][1] = .0;
#		} else {
#			$Mincoh[$ix][$iy][2] = $Mincoh[$ix][$iy][1] + ($int - $Mincoh[$ix][$iy][1]) / ($j+1);
#			$Sincoh[$ix][$iy][2] = $Sincoh[$ix][$iy][1] + ($int - $Mincoh[$ix][$iy][1]) * ($int - $Mincoh[$ix][$iy][2]);
#			$Mincoh[$ix][$iy][1] = $Mincoh[$ix][$iy][2];
#			$Sincoh[$ix][$iy][1] = $Sincoh[$ix][$iy][2];
#		}
#		
#		# coherent sum
#		$amp[1] = $inp[$index][3];
#		$amp[2] = $inp[$index][4];
#		$sumamp[$ix][$iy][1] += $amp[1];
#		$sumamp[$ix][$iy][2] += $amp[2];
#		
#		# compute running stddev according to Knuth's algorithm (https://www.johndcook.com/blog/standard_deviation/)
#		if ($j == 0) {
#			$Mcoh[$ix][$iy][1][1] = $amp[1];
#			$Scoh[$ix][$iy][1][1] = .0;
#			$Mcoh[$ix][$iy][2][1] = $amp[2];
#			$Scoh[$ix][$iy][2][1] = .0;
#		} else {
#			$Mcoh[$ix][$iy][1][2] = $Mcoh[$ix][$iy][1][1] + ($amp[1] - $Mcoh[$ix][$iy][1][1]) / ($j+1);
#			$Mcoh[$ix][$iy][2][2] = $Mcoh[$ix][$iy][2][1] + ($amp[2] - $Mcoh[$ix][$iy][2][1]) / ($j+1);
#
#			$Scoh[$ix][$iy][1][2] = $Scoh[$ix][$iy][1][1] + ($amp[1] - $Mcoh[$ix][$iy][1][1]) * ($amp[1] - $Mcoh[$ix][$iy][1][2]);
#			$Scoh[$ix][$iy][2][2] = $Scoh[$ix][$iy][2][1] + ($amp[2] - $Mcoh[$ix][$iy][2][1]) * ($amp[2] - $Mcoh[$ix][$iy][2][2]);
#			
#			$Mcoh[$ix][$iy][1][1] = $Mcoh[$ix][$iy][1][2];
#			$Scoh[$ix][$iy][2][1] = $Scoh[$ix][$iy][2][2];
#		}
#	}
