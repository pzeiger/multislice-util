#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;
use Data::Dumper;

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
	
	# build array of inputfiles
	if ( $setin == 1 ) {
		push @filein, $ARGV[$i];
		next;
	}
}

print "$lx, $ly, $accv, $nx, $ny, $detcangx, $detcangy, $innerdetang, $outerdetang, $setfort33\n";

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


my @header = ();
my @count = ();
my $ibeamymax = 0;

for my $j ( 0 .. $#filein ) {
    
    my @inp = ();
    my $iymin;
    my $ixmin;
    my $oldx;
    my $dat = [];
    
    print "opening input $filein[$j] \n";
    open my $in, '<', $filein[$j]
        or do{
	    warn "Could not open file $filein[$j] $!";
            next;
        };
    
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
	                if (! defined $ixmin or $inp[0] < $ixmin ) {
                        $ixmin = $inp[0];
                    }
	                if ( ! defined $iymin or $inp[1] < $iymin ) {
                        $iymin = $inp[1];
                    }
#                    print "$thetax $detcangx $innerdetang \n";
#                    print "$thetay $detcangy $innerdetang \n";
#                    print "quality of small angle approx ", abs($thetax-$thetax2), "\n";
                    my $thetadet = sqrt(($thetax-$detcangx)**2 + ($thetay-$detcangy)**2);
#                    print "$thetadet\n";
                    my @tmparr = ($inp[0], $inp[1], $inp[2], $inp[3], $thetadet);
                    push (@$dat, \@tmparr); 
                }
            }
        }
    }
    
    my $fileout = $filein[$j] . "_red";
    print "opening output $fileout\n";
    open my $out, '>', $fileout
        or do{
	    warn "Could not open file $fileout $!";
            next;
        };
    
    if ( $setfort33 ) {
        printf $out "# Reduced diffraction pattern\n"; 
        printf $out "# 1 -> px xpos\n"; 
        printf $out "# 2 -> px ypos\n"; 
        printf $out "# 3 -> real part of wave fcn\n";
        printf $out "# 4 -> imag part of wave fcn\n";
    }
    
    for my $row (@$dat) {
#        print "$row\n";
#        print Dumper(\@{$row});
        if ( defined $oldx and $row->[0] > $oldx ) {
             printf $out "\n";
        }
        if ( ! defined $oldx or $row->[0] > $oldx ) {
            $oldx = $row->[0];
        }
        if ( $setfort33 ) {
            printf $out "%4i %4i %+1.12e %+1.12e\n", ($row->[0]-$ixmin+1), ($row->[1]-$iymin+1), $row->[2], $row->[3];
        }
    }
    printf $out "\n";
    close($out);
}


