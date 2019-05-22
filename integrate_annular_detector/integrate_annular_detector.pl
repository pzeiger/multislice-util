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
my ($detcangx, $detcangy) = (0.0, 0.0);
my ($nbeamposx, $nbeamposy) = (0.0, 0.0);
my ($ly, $lx) = (1.0, 1.0); # in angstroms
my ($nx, $ny) = (0, 0);
my $accv = 100e3; # in V
my $setfort33 = 0;
my $setroiout = 0;
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
    # choose if to run in fort33 mode
	if ( $ARGV[$i] =~ /^-mode$/ ) {
		if ( $ARGV[$i+1] =~ /^fort33$/) {
		    $setfort33 = 1;
		}
		next;
	}
	# turn output of roi on
	if ( $ARGV[$i] =~ /^-roiout$/ ) {
        print "Output of ROI on\n";
		$setroiout = 1;
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
my $t = $accv * ELEMENTARY_CHARGE;                                  # kinetic energy, accv in keV
my $v = SPEED_OF_LIGHT * sqrt(1-(($r*$r)/(($r+$t)*($r+$t))));       # relativistic velocity
print $v / SPEED_OF_LIGHT;
my $mr = MASS_ELECTRON / sqrt(1-(($v*$v)/($csq)));                  # relativistic mass
my $lamb = 1e10 * PLANCK_CONSTANT / ($mr*$v);                       # wave length in angstroms
print "\n$lamb\n";
my $k = PI / $lamb;                                                 # length of k-vector
my $dkx = PI / ($lx);                                               # kx-spacing in 1/angstrom
my $dky = PI / ($ly);                                               # ky-spacing in 1/angstrom

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


if ( $setfort33 ) {
    printf $intout "# Integrated STEM intensity\n"; 
    printf $intout "# 1 -> incohint\n"; 
} else {
    printf $intout "# Integrated STEM intensity\n"; 
    printf $intout "# 1 -> x beampos\n";
    printf $intout "# 2 -> y beampos\n";
    printf $intout "# 3 -> incohint\n";
    printf $intout "# 4 -> uincohint\n";
    printf $intout "# 5 -> cohint\n";
    printf $intout "# 6 -> ucohint\n";
}


my @sumincoh = ();
my @usumincoh = ();
my @sumcoh = ();
my @usumcoh = ();
my @sumtds = ();
my @usumtds = ();
my @header = ();
my @count = ();

for my $j ( 0 .. $#filein ) {
    
    my @inp = ();
    my $ixmax = 0;
    
    print "opening $filein[$j] \n";
    open my $in, '<', $filein[$j]
        or do{
	    warn "Could not open file $filein[$j] $!";
            next;
        };
    
    my $beamposx = "";
    my $beamposy = "";
    if ( not $setfort33 ) {
        my @tmp = split('/', $filein[$j]);
        for my $i (0 .. $#tmp) {
            if ( $tmp[$i] =~ /^beampos_/ ) {
                my @tmp2 = split(/_x|_y/, $tmp[$i]);
                $beamposx = int($tmp2[1]);
                $beamposy = int($tmp2[2]);
                last;
            }
        } 
        if ( $beamposx eq "" and $beamposy eq "") {
            print "Error: Could not determine beam position. Please check file \"$filein[$j]\"!\n";
            die;
        }
        print "Current beam position: ($beamposx, $beamposy)\n";
    } else {
        $beamposx = 0;
        $beamposy = 0;
    }
    my $roiout;
    if ( $setroiout ) {
        my $roifileout = $filein[$j] . '_detc' . $detcangx . '-' . $detcangy . '_cang' . $innerdetang . '-'. $outerdetang . '_roi';
        
        open $roiout, '>', $roifileout
            or do{
                warn "Could not open file ${roifileout} $!\n";
                next;
            };
    }
    
    # initialize arrays for further use
    $count[$beamposx][$beamposy] = 0;
    $sumincoh[$beamposx][$beamposy] = 0.0;
    $usumincoh[$beamposx][$beamposy] = 0.0;
    $sumcoh[$beamposx][$beamposy] = 0.0;
    $usumcoh[$beamposx][$beamposy] = 0.0;
    $sumtds[$beamposx][$beamposy] = 0.0;
    $usumtds[$beamposx][$beamposy] = 0.0;
    
    while ( <$in> ) {
        if ( $_ =~ /^\s$/ ) {
            next;
        } elsif ( $_ =~ /^\#/) {
            push @header, $_;
            next;
        } else {
            @inp = grep { /\S/ } split(/\s+/, $_);
            
            my $thetax = ($inp[0]-$dpcnx)*$dkx/$k*1e3;  # in mrad
            my $thetax2 = asin(($inp[0]-$dpcnx)*$dkx/$k)*1e3;  # in mrad
#            print $thetax, $detcangx, $innerdetang, "\n";
            if ( abs($thetax-$detcangx) <= $outerdetang ) {
                my $thetay = ($inp[1]-$dpcny)*$dky/$k*1e3;  # in mrad
                if ( abs($thetay-$detcangy) <= $outerdetang ) {
#                    print "thetax: ", $thetax, " ", $inp[0], "\n";
#                    print "thetay: ", $thetay, " ", $inp[1], "\n";
                    if ( $inp[0] > $ixmax and $ixmax ne 0 and $setroiout ) {
                         printf $roiout "\n";
                    }
	                if ( $inp[0] > $ixmax ) {
                        $ixmax = $inp[0];
                    }
#                    print "$thetax $detcangx $innerdetang \n";
#                    print "$thetay $detcangy $innerdetang \n";
#                    print "quality of small angle approx ", abs($thetax-$thetax2), "\n";
                    my $thetadet = sqrt(($thetax-$detcangx)**2 + ($thetay-$detcangy)**2);
#                    print "$thetadet\n";
                     
                    if ( $setfort33 ) {
                        if ( $thetadet >= $innerdetang and $thetadet <= $outerdetang ) {
#                            print "$count ok \n";
                            $count[$beamposx][$beamposy] += 1;
                            $sumincoh[$beamposx][$beamposy] += $inp[2]*$inp[2] + $inp[3]*$inp[3];
                            if ( $setroiout ) {
                                printf $roiout "%4i %4i %+1.12e %+1.12e\n", $inp[0], $inp[1], $inp[2], $inp[3];
                            }
                        } else {
                            if ( $setroiout ) {
                                printf $roiout "%4i %4i %+1.12e %+1.12e\n", $inp[0], $inp[1], 0.0, 0.0;
                            }
                        }
                    } else {
                        if ( $thetadet >= $innerdetang and $thetadet <= $outerdetang ) {
#                            print "$count ok \n";
                            $count[$beamposx][$beamposy] += 1;
                            $sumincoh[$beamposx][$beamposy] += $inp[2];
                            $usumincoh[$beamposx][$beamposy] += $inp[3]*$inp[3];
                            $sumcoh[$beamposx][$beamposy] += $inp[4];
                            $usumcoh[$beamposx][$beamposy] += $inp[5]*$inp[5];

                            $sumtds[$beamposx][$beamposy] += ($inp[2] - $inp[4]);
                            $usumtds[$beamposx][$beamposy] += ($inp[3]**2 + $inp[5]**2);
                            
                            if ( $setroiout ) {
                                print "$setroiout\n";
                                printf $roiout "%4i %4i %1.12e %1.12e %1.12e %1.12e\n", $inp[0], $inp[1], $inp[2], $inp[3], $inp[4], $inp[5];
                            }
                        } else {
                            if ( $setroiout ) {
                                printf $roiout "%4i %4i %1.12e %1.12e %1.12e %1.12e\n", $inp[0], $inp[1], 0.0, 0.0, 0.0, 0.0;
                            }
                        }
                    }
                }
            }
        }
    }
    
    if ( $setroiout ) {
        printf $roiout "\n";
        close($roiout);
    }
    
    if ( $setfort33 ) {
        printf $intout "%1.12e\n", $sumincoh[$beamposx][$beamposy];
    }
}

if ( not $setfort33 ) {
    for my $ibx (0 .. $#sumincoh) {
        for my $iby (0 .. $#{$sumincoh[$ibx]}) {
            if ( defined($sumincoh[$ibx][$iby]) ) {
                printf $intout "%2i %2i %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e\n", $ibx, $iby, $sumincoh[$ibx][$iby], sqrt($usumincoh[$ibx][$iby]), $sumcoh[$ibx][$iby], sqrt($usumcoh[$ibx][$iby]), $sumtds[$ibx][$iby], sqrt($usumtds[$ibx][$iby]);
            } else {
                print "Warning: data missing for beam position ($ibx, $iby). Printing zeros...\n";
                printf $intout "%2i %2i %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e %1.12e\n", $ibx, $iby, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
            }
        }
        printf $intout "\n";
    }
}


