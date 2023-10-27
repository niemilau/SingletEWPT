#!/bin/bash
# This script is file structure dependent. test/mh2_$mh2/0~#
# This script is used to collect the transitions in mh2_$mh2, where $mh2 is in the range of [mh2_min, mh2_max].
# combine these transition results in the current directory file named results.
if [ "$#" -ne 2 ]; then
    echo "Usage: <script name> mh2_min mh2_max"
	exit 1
fi

mh2_min=$1
mh2_max=$2

transitions="transitions.dat"
results="results"

if [[ ! -f "$results" ]]; then
    touch $results
fi

pathResults="$(pwd)/$results"
echo $pathResults

for dir in */
do
    if [[ $dir != mh2-* ]]; then
        continue
    fi

    mh2=${dir#mh2-}
    mh2=${mh2%/*}
    
    if (( $(echo "$mh2 > $mh2_max" |bc -l) || $(echo "$mh2 < $mh2_min" |bc -l) )); then
		continue
	fi

    cd $dir
    

    for subdir in */
    do
        cd $subdir


        if [[ ! -f "$transitions" ]]; then
      		echo "Skipping directory $dir$subdir"
            cd ..
      		continue
    	fi

        cat "$transitions" >> "$pathResults"
        cd ..
    done
    cd ..
done

