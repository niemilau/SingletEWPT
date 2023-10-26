#!/bin/bash
# The directory structure should be like run/mh2-$mh2/<splitted a2>, which is based on the PrepareBigScan.sh
# In the directory run/, call this script and it will go through all leaf directories.
# If the leaf directory has submit_scanner.job and parameters file for SingletEWPT, submit the job in the INPAC-Cluster
scanner_job='submit_scanner.job'
paramfile='parameters'
for dir in */
do
    if [[ $dir == "*/" ]]; then
        continue
    fi
    cd $dir

    for subdir in */
    do
        if [[ $subdir == "*/" ]]; then
            continue
        fi

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