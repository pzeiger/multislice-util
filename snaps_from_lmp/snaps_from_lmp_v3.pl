#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my %atomlist;
my $atomstyle;
my $centroidfile;
my $centroidin;
my $formatin = "MD";
my $in;
my $lengthunit = "angstrom";
my $mode = "clmd";
my $nhead;
my ($setin, $setatom) = (0, 0,);
my $setfracin = 0;
my $setfracout = 1;
my $setsymout = 0;
my $trlcorr = 0;
my $atomtype;
my @atomtypes;

#my $snapcount = 1;
#my $outevery = 1;
# start output after timesteps
my $nfirst = 0;
# -1 -> sample until EOF, other wise give end time in time units
my $nlast = -1;
# sample every
my $nsample = 1;

# Loop over arguments and place them into arrays and vars for further use
foreach my $i ( 0 .. $#ARGV ) {
    
    # Specify LAMMPS atom type. Currently supported: "full" and "atomic"
    if ( $ARGV[$i] =~ /^-atomstyle$/ ) {
        ($setin, $setatom) = (0, 0,);
        $atomstyle = $ARGV[$i+1];
        if ( $atomstyle ne "full" and $atomstyle ne "atomic" ) {
            print "atom style $atomstyle not recognized.\n";
            die;
        }
        next;
    }
    
    # Recognize which input file is to be used
    if ( $ARGV[$i] =~ /^-in$/ ) {
        ($setin, $setatom) = (1, 0,);
        next;
    }
    
    # Set sampling properties
    if ( $ARGV[$i] =~ /^-sample$/ ) {
	    $nfirst = $ARGV[$i+1];
	    $nsample = $ARGV[$i+2];
	    $nlast = $ARGV[$i+3];
	    if ($nfirst and ($nfirst < $nlast or $nlast < 0) and $nsample > 0) {
            next;
	    } else {
            print('False input for -sample');
            die;
        }
    }
    
    # Match atom type from lammps with nuclear charge and atom style
    if ( $ARGV[$i] =~ /^-matchatom$/ ) {
        ($setin, $setatom) = (0, 1,);
	if (defined $atomtype) {
            $atomtype += 1;
        } else {
            $atomtype = 0;
	}
        next;
    }
    
    # Specify the unit of lenght in the input file. Currently supported are "angstrom" and "bohr"
    if ( $ARGV[$i] =~ /^-lengthunit$/ ) {
        ($setin, $setatom) = (0, 0,);
        $lengthunit = $ARGV[$i+1];
        if ($lengthunit eq "angstrom" or $lengthunit eq "bohr") {
            print "Unit of length is set to $lengthunit\n";
        } else {
            print "Unit of length not recognized\n";
            die;
        }
        next;
    }
    
    # input centroid file for COM correction in output of PIMD simulation
    if ( $ARGV[$i] =~ /^-fcentroid$/ ) {
        ($setin, $setatom) = (0, 0,);
        $centroidfile = $ARGV[$i+1];
        next;
    }
    
    # Specify operation mode. Currently supported: "pimd" and "clmd"
    if ( $ARGV[$i] =~ /^-mode$/ ) {
        ($setin, $setatom) = (0, 0,);
        $mode = $ARGV[$i+1];
        next;
    }
    
    # specify that input file uses fractional coordinates
    if ( $ARGV[$i] =~ /^-fracin$/ ) {
        ($setin, $setatom) = (0, 0,);
        $setfracin = 1;
        next;
    }

    # specify that input file uses fractional coordinates
    if ( $ARGV[$i] =~ /^-nofracout$/ ) {
        ($setin, $setatom) = (0, 0,);
        $setfracout = 0;
        next;
    }
    
    # specify that simulation box in output file shall be symmetric around origin
    if ( $ARGV[$i] =~ /^-symout$/ ) {
        ($setin, $setatom) = (0, 0,);
        $setsymout = 1;
        next;
    }
    
    # Specify the input format. Currently supported are "MD", "lmp" and "xyz"
    # Variable $nhead specifies the numebr of lines containing 
    # environment information (non-position information) per time step
    # in lammps output
    if ( $ARGV[$i] =~ /^-format$/ ) {
        ($setin, $setatom) = (0, 0,);
        $formatin = $ARGV[$i+1];
        if ($formatin eq "lmptrj") {
            $nhead = 8;
            print "Input file format is lmptrj\n";
        } elsif ($formatin eq "lmpdata") {
            $nhead = 2;
            print "Input file format is lmpdata\n";
        } elsif ($formatin eq "ipixyz") {
            $nhead = 2;
            print "Input file format is ipixyz\n";
        } else {
            print "Input file format not recognized\n";
            die;
        }
        next;
    }
    
    # Turn the center of mass (com) translation correction on (1) or off (all other values)
    if ( $ARGV[$i] =~ /^-trlcorr$/ ) {
        ($setin, $setatom) = (0, 0,);
        print "COM translation correction turned on\n";
        $trlcorr = 1;
        next;
    }
    
    # set all the calculation parameters
    if ( $setin == 1 ) {
        open $in, '<', $ARGV[$i];
        $setin = 0;
        next;
    }
    if ( $setatom == 1 ) {
#	print "\n$ARGV[$i]\n";
	if ($ARGV[$i] =~ /:\z/ ) {
            $atomtypes[$atomtype] = int(substr($ARGV[$i], 0, length($ARGV[$i])-1));
        } else {
            $atomlist{$ARGV[$i]} = $atomtypes[$atomtype];
#	    print "$atomlist{$ARGV[$i]}\n";
	}
#	print "$atomtypes[$atomtype]\n";
#	print "$atomtype\n\n";
        next;
    }
}

#print Dumper(\%atomlist);
#die;

# Print some info and open pimd centroid file
if ( not $trlcorr) {
    print "COM translation correction turned off\n";
}

if ( $mode =~ /pimd/ ) {
    print "Operating in PIMD mode\n";
    open $centroidin, '<', $centroidfile;
} else {
    print "Operating in CLMD mode\n";
}

# create initial structure array; if given use first snapshot of the
# input file if not extra file ($ARGV[2]) is specified 

my $line = <$in>;
my ($tstep0, $natom0, $snap0, $dim0, $ixyz0, $typecol0) = readTimestepData($in,$line,$nhead,$formatin,$lengthunit,$setfracin,$atomstyle);

my @snap0 = @$snap0;
my @dim0 = @$dim0;
my @ixyz0 = @$ixyz0;


my @centroid = ();
my $centroid;

if ( $mode =~ /pimd/ && $trlcorr == 1 ) {
    my $line = <$centroidin>;
    
    (my $tstep, my $natom, $centroid, my $dim, my $ixyz) = readTimestepDataMD($centroidin,$line,$nhead,$formatin,$lengthunit,$setfracin,$atomstyle);
    
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

print "\nTIMESTEP: $tstep0\n";
print "CENTRE OF MASS: @com0\n";

if ( $mode !~ /pimd/ ) {
    my $fileout = "snapshot" . $tstep0;
    open my $out,'>', $fileout;
    
    printf $out "%f %f %f\n", @dim0;
    printf $out "%i F\n", $natom0;
    foreach ( @snap0 ) {
        
        if ( @{$_}[$ixyz0[0]] < 0 ) {@{$_}[$ixyz0[0]] += 1;}
        if ( @{$_}[$ixyz0[1]] < 0 ) {@{$_}[$ixyz0[1]] += 1;}
        if ( @{$_}[$ixyz0[2]] < 0 ) {@{$_}[$ixyz0[2]] += 1;}
        
        if ( @{$_}[$ixyz0[0]] > 1 ) {@{$_}[$ixyz0[0]] -= 1;}
        if ( @{$_}[$ixyz0[1]] > 1 ) {@{$_}[$ixyz0[1]] -= 1;}
        if ( @{$_}[$ixyz0[2]] > 1 ) {@{$_}[$ixyz0[2]] -= 1;}
        
	if ( $setsymout ) {
            @{$_}[$ixyz0[0]] -= 0.5;
            @{$_}[$ixyz0[1]] -= 0.5;
            @{$_}[$ixyz0[2]] -= 0.5;
	}
        
	if ( not $setfracout ) {
            @{$_}[$ixyz0[0]] *= $dim0[0];
            @{$_}[$ixyz0[1]] *= $dim0[1];
            @{$_}[$ixyz0[2]] *= $dim0[2];
	}
        
        my $atomtype = @{$_}[$typecol0];
        printf $out "%2i %3.16e %3.16e %3.16e\n", $atomlist{$atomtype}, @{$_}[$ixyz0[0]], @{$_}[$ixyz0[1]], @{$_}[$ixyz0[2]];
    }
}
else {
    my $fileout = "snapshot" . $tstep0;
    open my $out,'>', $fileout;
    
    printf $out "%f %f %f\n", @dim0;
    printf $out "%i F\n", $natom0;
    foreach ( @centroid ) {
        my $atomtype = @{$_}[$typecol0];
        printf $out "$atomlist{$atomtype} @{$_}[$ixyz0[0]] @{$_}[$ixyz0[1]] @{$_}[$ixyz0[2]] \n";
    }
}


# Print all other snapshots to multislice input files
while ( my $line = <$in> ) {
    
    my ($tstep, $natom, $snap, $dim, $ixyz, $typecol) = readTimestepData($in,$line,$nhead,$formatin,$lengthunit,$setfracin,$atomstyle);
    
    if (not $tstep) {
        last;
    } elsif ($tstep < $nfirst) {
		print "time step not sampled\n";
        next;
    } elsif ($tstep > $nlast and $nlast > 0) {
		print "time step not sampled\n";
        last;
    } elsif ($tstep % $nsample != 0) {
		print "time step not sampled\n";
        next;
    }
    
    my @snap = @$snap;
    my @dim = @$dim;
    
    my @ixyz = @$ixyz;
    
    my @comcorr = ();
    
    # correct for atoms reentering simulation box on the opposite site
    # due to periodic boundary conditions
    foreach my $i ( 0 .. $#snap ) {
        foreach my $j ( $ixyz[0] .. $ixyz[2] ) {
#            print "@ixyz\n";
#            print "$i $j $snap[$i][$j]-$snap0[$i][$j]\n";
            my $disp = $snap[$i][$j] - $snap0[$i][$j];
#            print "$i $j =  $disp\n";
            
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
        (my $tstep, my $natom, $centroid, my $dim, my $ixyz) = readTimestepData($centroidin,$line,$nhead,$formatin,$lengthunit,$setfracin,$atomstyle);
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
        
        if ( @{$_}[$ixyz[0]] < 0 ) {@{$_}[$ixyz[0]] += 1;}
        if ( @{$_}[$ixyz[1]] < 0 ) {@{$_}[$ixyz[1]] += 1;}
        if ( @{$_}[$ixyz[2]] < 0 ) {@{$_}[$ixyz[2]] += 1;}
        
        if ( @{$_}[$ixyz[0]] > 1 ) {@{$_}[$ixyz[0]] -= 1;}
        if ( @{$_}[$ixyz[1]] > 1 ) {@{$_}[$ixyz[1]] -= 1;}
        if ( @{$_}[$ixyz[2]] > 1 ) {@{$_}[$ixyz[2]] -= 1;}
        
        
        if ( $setsymout ) {
            @{$_}[$ixyz[0]] -= 0.5;
            @{$_}[$ixyz[1]] -= 0.5;
            @{$_}[$ixyz[2]] -= 0.5;
        }
        
        if ( not $setfracout ) {
            @{$_}[$ixyz[0]] *= $dim[0];
            @{$_}[$ixyz[1]] *= $dim[1];
            @{$_}[$ixyz[2]] *= $dim[2];
        }
        
        my $atomtype = @{$_}[$typecol];    
        printf $out "$atomlist{$atomtype} @{$_}[$ixyz[0]] @{$_}[$ixyz[1]] @{$_}[$ixyz[2]]\n";
    }
}




# subroutine that reads data for 1 time step from a LAMMPS dump file formatted as:
# 
# 1. MD file
#    ------BEGIN FILE------
#    ITEM: TIMESTEP
#    0
#    ITEM: NUMBER OF ATOMS
#    34680
#    0.0  62.426 xlo xhi
#    0.0  61.271 ylo yhi
#    0.0 300.165 zlo zhi
#    ITEM: ATOMS id type xs ys zs
#    1 1 0.0000000 0.0000000 0.0000000
#    2 1 0.0333333 0.0000000 0.0098039
#    ...
#    ------END FILE ------
#
#
# 2. xyz file
#    ------BEGIN FILE------
#    34680
#    # CELL(abcABC):  117.96804   115.78541   567.22963    90.00000    90.00000    90.00000  cell{atomic_unit}  Traj: positions{atomic_unit} Step:           0  Bead:       0 
#    Gd  0.00000e+00  0.00000e+00  0.00000e+00
#    Gd  3.93252e+00  0.00000e+00  5.56146e+00
#    ...
#    ------END FILE ------
#   
#   
# 3. lmp lammps datafile
#    ------BEGIN FILE------
#    Multilayered hBN 12x21x45
#   
#    45360 atoms
#    2 atom types
#    
#    -52.200000 52.200000 xlo xhi
#    -52.740947 52.740947 ylo yhi
#    -149.850000 149.850000 zlo zhi
#    
#    Masses
#    1 10.810000
#    2 14.007000
#    
#    Atoms
#    1 1 1 0.000000 -26.100000 -26.370474 -73.260000
#    2 1 2 0.000000 -26.100000 -25.951895 -73.260000
#    ...
#   ------END FILE ------
#
#    
# takes 4 arguments. Current input line $line, input file identifier $in, number
# of lines $nhead composing the header, and the input file format (currently .xyz, .MD) --> see examples


sub readTimestepData {
    my $in = $_[0];
    my $line = $_[1];
    my $nhead = $_[2];
    my $formatin = $_[3];
    my $lengthunit = $_[4];
    my $fractional = $_[5];
    my $atomstyle = $_[6];
    
    
    # arrays for unit cell dimensions, atom data and temporary auxiliary purposes
    my ($tstep, $natom);
    my @dim = ();
    my @snap = ();
    my %data;
    my %head;
    my $idcol;
    my $typecol;
	my $cellunit = $lengthunit;
	my $posunit = $lengthunit;
    
    my @ixyz;
    
    
    # Read header
    chomp($line);
    if ($line ne '' ) {
        $head{0} = $line;
        
        if ( $formatin =~ /lmptrj/ ) {
            foreach my $i (1 .. $nhead) {
                $line = <$in>;
                chomp($line);
                $head{$i} = $line;
#            print "$i $head{$i} \n";
            }
            
            # 1st line of the header contains time step info;
            # 3rd line the number of atoms
            $tstep = (split ' ', $head{1})[0];
            $natom = (split ' ', $head{3})[0];
            
            my @tmp = split ' ', $head{$nhead};
            
            # determine column number of id
            foreach my $i (0 .. $#tmp) {
                if ( $tmp[$i] =~ /id/ ) {
                    $idcol = $i;
                }
                if ( $tmp[$i] =~ /type/ ) {
                    $typecol = $i;
                }
            }
            $typecol = $typecol - $idcol;
            
            foreach my $i (0 .. $#tmp ) {
                if ( $tmp[$i] =~ /xs/ ) {
                    $ixyz[0] = $i - $idcol;
                }
                if ( $tmp[$i] =~ /^x$/ ) {
                    $ixyz[0] = $i - $idcol;
                }
                if ( $tmp[$i] =~ /ys/ ) {
                    $ixyz[1] = $i - $idcol;
                }
                if ( $tmp[$i] =~ /^y$/ ) {
                    $ixyz[1] = $i - $idcol;
                }
                if ( $tmp[$i] =~ /zs/ ) {
                    $ixyz[2] = $i - $idcol;
                }
                if ( $tmp[$i] =~ /^z$/ ) {
                    $ixyz[2] = $i - $idcol;
                }
            }
            
            foreach my $i ( ($nhead-3) .. ($nhead-1) ) {
                my @tmp = split ' ', $head{$i};
                push @dim, ( $tmp[1] - $tmp[0] );
            }
            
            foreach ( 1 .. $natom ) {
                $line = <$in>;
                chomp($line);
                @tmp = split ' ', $line;
                $data{$tmp[0]} = $line;
                push @snap, [ @tmp ];
            }
            
        } elsif ( $formatin =~ /xyz/ ) {
            
            foreach my $i (1 .. ($nhead-1)) {
                $line = <$in>;
                chomp($line);
                $head{$i} = $line;
#            print "$i $head{$i} \n";
            }
            
            $natom = (split ' ', $head{0})[0];
            my @tmp = split ' ', $head{1};
            @dim = @tmp[2,3,4];
            $tstep = $tmp[9];
            
			if ($tmp[12] =~ /position{atomic_unit}/) {
				$cellunit = "bohr";
			}
			if ($tmp[13] =~ /cell{atomic_unit}/) {
				$cellunit = "bohr";
			}
            foreach (1 .. $natom) {
                $line = <$in>;
                chomp($line);
                push @snap, [split ' ', $line];
            }
            
            # 2nd, 3rd, and 4th column contain the x-, y- and z-positions
			$typecol = 0;
            @ixyz = (1, 2, 3);
            
        } elsif ( $formatin =~ /lmpdata/ ) {
            my $i = 0;
            while ($line = <$in>) {
                chomp($line);
                print "$line\n";
                my @tmp = split ' ', $line;
                if ( !@tmp ) {
                    next;
                }
                if ( $tmp[0] =~ /[Aa]toms/ ) {
                    last;
                } elsif ($tmp[0] =~ /[Mm]asses/ ) {
                    foreach (1 .. $head{'Ntypes'}) {
                        $line = <$in>;
                        chomp($line);
                        if ($line =~ //) {
                    	    $line = <$in>;
                            chomp($line);
                        }
                        my @tmp2 = split ' ', $line;
                        print @tmp2;
                        $head{"mass$tmp2[0]"} = $tmp2[1];
                    }
                } elsif ($tmp[1] =~ /[Aa]toms/ ) {
                    $head{'Natoms'} = $tmp[0];
                } elsif ($tmp[1] =~ /[Aa]tom/ and $tmp[2] =~ /[Tt]ypes/) {
                    $head{'Ntypes'} = $tmp[0];
                } elsif ($tmp[2] =~ /xlo/ and $tmp[3] =~ /xhi/ ) {
                	print "$tmp[0]\n";
                	print "$tmp[2]\n";
                    $head{'xlo'} = $tmp[0];
                    $head{'xhi'} = $tmp[1];
            	print "$head{'xlo'}\n";
                } elsif ($tmp[2] =~ /ylo/ and $tmp[3] =~ /yhi/ ) {
                    $head{'ylo'} = $tmp[0];
                    $head{'yhi'} = $tmp[1];
            	print "$head{'xlo'}\n";
                } elsif ($tmp[2] =~ /zlo/ and $tmp[3] =~ /zhi/ ) {
                    $head{'zlo'} = $tmp[0];
                    $head{'zhi'} = $tmp[1];
            	print "$head{'xlo'}\n";
                }
                $i += 1;
            }
            print Dumper(\%head);
#            print "$head{'xlo'}\n";
            $natom = $head{'Natoms'};
#            print $head{'xhi'};
            @dim = ( $head{'xhi'}-$head{'xlo'}, $head{'yhi'}-$head{'ylo'}, $head{'zhi'}-$head{'zlo'} );
            $tstep = 0;
            
            foreach (1 .. $natom) {
                $line = <$in>;
                chomp($line);
                
                if ($line =~ //) {
                    $line = <$in>;
                    chomp($line);
                }
                push @snap, [split ' ', $line];
            }
            
            # 2nd, 3rd, and 4th column contain the x-, y- and z-positions
            if ($atomstyle =~ /full/) {
                $typecol = 2;
                @ixyz = (4, 5, 6);
            } else {
				print("Unsupported atomtype");
				die;
			}
            
            for my $i (0 .. $natom-1) {
                $snap[$i][$ixyz[0]] -= $head{'xlo'};
                $snap[$i][$ixyz[1]] -= $head{'ylo'};
                $snap[$i][$ixyz[2]] -= $head{'zlo'};
            }
        }
        
        # Convert cell dimension to Angstroms
        if ($cellunit =~ /bohr/) {
            print "Unit of length for cell dimensions is $cellunit... Converting to angstrom\n";
			# numerical value of bohr taken from ipi
            $dim[0] *= 0.529177249;
            $dim[1] *= 0.529177249;
            $dim[2] *= 0.529177249;
        } elsif ($cellunit ne "angstrom") {
            print "Unsupported unit. Terminating....\n";
	        die;
        }

        if ($posunit =~ /bohr/) {
            print "Unit of length for atomic positions is $posunit... Converting to angstrom\n";
			# numerical value of bohr taken from ipi
            for my $i (0 .. $natom-1) {
                $snap[$i][$ixyz[0]] *= 0.529177249;
                $snap[$i][$ixyz[1]] *= 0.529177249;
                $snap[$i][$ixyz[2]] *= 0.529177249;
			}
        } elsif ($posunit ne "angstrom") {
            print "Unsupported unit. Terminating....\n";
	        die;
        }
		
        # convert to fractional coordinates if necessary
        if (not $fractional) {
            for my $i (0 .. $natom-1) {
                $snap[$i][$ixyz[0]] /= $dim[0];
                $snap[$i][$ixyz[1]] /= $dim[1];
                $snap[$i][$ixyz[2]] /= $dim[2];
            }
        }
    }
    
    print "@dim \n$tstep\n";
    # \@ creates a reference to the array -> see:
    # http://perlmeme.org/faqs/perl_thinking/returning.html	
    return($tstep, $natom, \@snap, \@dim, \@ixyz, $typecol);
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

