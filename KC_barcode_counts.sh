#!/bin/bash

# First argument is input file (fastq), second argument is output file root,
# third argument is length of labeling BC (usually 9 or 12)
ARGS=("$@")

# Check for correct number of arguments
if [ "$#" -ne 3 ]; then
	echo "Usage: [input file] [output file root] [BC length]"
	exit 1
fi

BCLENGTH=${ARGS[2]}

# Pulls lines from fastq file that contain sequence and puts in a temporary file, changed length to 101 from 50
TMPFILE="tmp"$(date +%s)
grep -e "[AGCTN]\{101,\}" ${ARGS[0]} > $TMPFILE

# Multiplex BCs used in the sample and numbers that refer to them
# FILL IN HERE

#MPBC=('GCTCGAT' 'TAGACTAT' 'CGCTACCCT' 'ATAGTGGACA' 'GTCAGTAGGTA')
MPBC=('GCTCGAT' 'TAGACTAT' 'CGCTACCCT' 'ATAGTGGACA' 'GTCAGTAGGTA')

BCNUM=(01 02 03 04 05)

# MIDDLE is sequence between multiplex BC and labeling BC
#MIDDLE='TCGAG' #Used for pCDNA5-RE_ARRAY-1
MIDDLE='TCGAG' #Used for pCDNA5-RE_ARRAY-1

# Loops through each multiplex BC used in the sample
for k in $(seq 0 $(expr ${#MPBC[*]} - 1))
do
	echo "Processing BC "${BCNUM[$k]}
	CURRENTFILE=${ARGS[1]}'_'${BCNUM[$k]}
	
	# Pulls the labeling BC sequence from each line and outputs it, only if 
	# the sequence has the proper middle sequence and current multiplex BC
	sed -n 's/.*'${MPBC[$k]}$MIDDLE'\([ACGT]\{'$BCLENGTH'\}\).*/\1/p' $TMPFILE > $CURRENTFILE
	
	# Sorts the labeling BCs (necessary for uniq command) and then collapses
	# them by the unique BC sequence with counts of how many times each
	# labeling BC appears in the file
	sort $CURRENTFILE | uniq -c >$CURRENTFILE"_counts"
done

# Removes temporary file of just full sequences
rm $TMPFILE


