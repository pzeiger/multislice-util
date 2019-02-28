#!/bin/bash

DIRIN=$1

PREPDIR="prep_snap"
AVGDIR="avg_snaps"
PDIR=$(pwd)

while read -r line
do
    WDIR[$i]=${line//[[:blank:]]/}
    i="$((i+1))"
done < "$DIRIN"

n=${#WDIR[@]}

for (( c=0; c<$n; c++ )); do
	cd ${WDIR[$c]}
	
	mkdir -p $AVGDIR
	cd $AVGDIR

	#(cat $PDIR/${WDIR[$c]}/run*/snaps_run_*.in | xargs -n 1 -I {} 
	nohup $PDIR/avg_diffpatt_over_snaps.pl -in $PDIR/${WDIR[$c]}/*_mult/fort.33 -outevery 10 > avg.out &
	
	
    cd $PDIR
done



