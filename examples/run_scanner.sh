#!/bin/bash
# The directory structure should be like run/mh2-$mh2/<splitted a2>, which is based on the PrepareBigScan.sh
# In the directory run/, call this script with three arguments mh2_min, mh2_max and mh2_delta, 
# it will go through all directories in the mh2-$mh2 with step mh2_delat from [mh2_min, mh2_max].
# If the leaf directory has submit_scanner.job and parameters file for SingletEWPT, submit the job in the INPAC-Cluster
if [ "$#" -ne 3 ]; then
    echo "Usage: <script name> mh2_min mh2_max mh2_delta"
	exit 1
fi

mh2_min=$1
mh2_max=$2
mh2_delta=$3
scanner_job='submit_scanner.job'
paramfile='parameters'


for mh2 in $(seq $mh2_min $mh2_delta $mh2_max);
do
    dir="mh2-$mh2"
    if [ ! -d "$dir" ]; then
  		echo "Error: No directory $dir"
        continue
    fi
    cd $dir
    for subdir in */
    do
        cd $subdir

        if [[ ! -f "$scanner_job" || ! -f "$paramfile" ]]; then
  		    echo "Skipping directory $dir$subdir"
            cd ..
  		    continue
	    fi
        condor_submit $scanner_job
        cd ..
    done
    cd ..
done








: <<'END_COMMENT'

# This is the version to rm -rf some files

if [ "$#" -ne 3 ]; then
    echo "Usage: <script name> mh2_min mh2_max mh2_delta"
	exit 1
fi

mh2_min=$1
mh2_max=$2
mh2_delta=$3
scanner_job='submit_scanner.job'
paramfile='parameters'


for mh2 in $(seq $mh2_min $mh2_delta $mh2_max);
do
    dir="mh2-$mh2"
    if [ ! -d "$dir" ]; then
  		echo "Error: No directory $dir"
        continue
    fi
    cd $dir
    for subdir in */
    do
        cd $subdir
        if [[ ! -f "$scanner_job" || ! -f "$paramfile" ]]; then
  		    echo "Skipping directory $dir$subdir"
            cd ..
  		    continue
	    fi
        rm -rf temperature_data.dat
        cd ..
    done
    cd ..
done
END_COMMENT