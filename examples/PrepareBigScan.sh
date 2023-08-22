#!/bin/bash

### Use this script for splitting a large scan into smaller pieces (that can then be ran in parallel on a cluster)
# NB: Needs a template 'parameters' file with correct options etc in the directory where you run this

if [ "$#" -ne 3 ]; then
    echo "Usage: <script name> mh2_min mh2_max mh2_delta"
	exit 1
fi

PARAMFILE=parameters

if [ ! -f "$PARAMFILE" ]; then
	echo "Error: no input file '$PARAMFILE'"
	exit
fi

## additional files to copy for each run (eg. job submit file for slurm)
additional_files=()
additional_files+=("submit_scanner.job")

# Range for mh2 (will make a separate folder for each)

mh2_min=$1
mh2_max=$2
mh2_delta=$3

# For other parameters: will either use the ranges given in PARAMFILE,
# or read them from files like 'range_a2' if they exist

## Further split each mh2 run? 
# For this cookiecutter script I will just split the a2 range, others are kept intact
# Will also require a custom a2 range to be given
subranges_a2=4 # 1 = no a2 splitting
RANGEFILE_a2="range_a2"
SUBFILE_DIR="subfiles_a2"

## Store subfile names in an array for use in the mh2 loop
subfile_names=()

if [[ ! -f "$RANGEFILE_a2" && $subranges_a2 -gt 1 ]]; then
	echo "Error: no input file '$RANGEFILE_a2', can't split"
	exit
else
	## prepare smaller range_a2 files
	lines_total=$(wc -l < $RANGEFILE_a2)
	lines_per_subfile=$(($lines_total / $subranges_a2))
	remainder=$(($lines_total % $subranges_a2))
	if [ $lines_per_subfile -lt 1 ]; then
		echo "Error: can't divide $lines_total a2 values to $subranges_a2 smaller sets"
		exit
	fi
	# clear old files from subdir
	rm -rf $SUBFILE_DIR
	mkdir -p $SUBFILE_DIR
	split -d -l "$lines_per_subfile" "$RANGEFILE_a2" "$SUBFILE_DIR/range"
	# this produces files range00, range01 etc
	## NB: If the division was uneven, split will produce 1 extra file with less entries. Are we happy with this??
	if [ $remainder -gt 0 ]; then
		subranges_a2=$subranges_a2+1
	else
		subranges_a2=$subranges_a2
	fi


	echo "==== a2 subranges ===="
	for ((i=0; i<$subranges_a2; i++)); do
    	subfile="$SUBFILE_DIR/range0$i"
		echo "File $i:"
		echo $(cat $subfile)
		subfile_names+=($subfile)
	done
fi



for mh2 in $(seq $mh2_min $mh2_delta $mh2_max);
do

	for ((i=0; i<$subranges_a2; i++)); do 

		DIRNAME="mh2-$mh2/$i"

		if [ ! -d "$DIRNAME" ]; then
		
			mkdir -p $DIRNAME
			
			## Change mass in the parameters file. I want to keep mh2 fixed in each run => set min=max etc
			
			# (\s in sed will match whitespaces including tab)
			sed -i "s/mh2_min\s.*/mh2_min \t$mh2/g" $PARAMFILE
			sed -i "s/mh2_max\s.*/mh2_max \t$mh2/g" $PARAMFILE
			sed -i "s/mh2_delta\s.*/mh2_delta \t1/g" $PARAMFILE
		
			## Copy range_a2 file
			cp "${subfile_names[$i]}" "$DIRNAME/range_a2"

			## Copy scanning ranges for other parameters if they were given
			[ -f "range_b3" ] && cp "range_b3" "$DIRNAME/"
			[ -f "range_b4" ] && cp "range_b4" "$DIRNAME/"
			[ -f "range_sinTheta" ] && cp "range_sinTheta" "$DIRNAME/"
			[ -f "range_T" ] && cp "range_T" "$DIRNAME/"

			cp "$PARAMFILE" "$DIRNAME/" 

			# copy additional files
			for fname in "${additional_files[@]}"; do

				## For slurm job files: change job name
				if [[ -f "$fname" && "$fname" == "submit_scanner.job" ]]; then
					sed -i "s/--job-name=.*/--job-name=m$mh2-$i/g" submit_scanner.job
				fi
				
				[ -f "$fname" ] && cp "$fname" "$DIRNAME/"
			done
		
		else
			echo 'Error: directory' $DIRNAME 'exists.'
		fi

	done ## a2
done ## mh2
