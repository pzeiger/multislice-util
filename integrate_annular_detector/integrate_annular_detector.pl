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


