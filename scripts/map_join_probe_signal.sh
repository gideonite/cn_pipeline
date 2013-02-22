#!/bin/bash

# runs the function join_probe_signal from probe_utils.py over all the files in a
# directory puts the output either in the specified directory or in the randomly
# generated directory </tmp/cbs_in.$RANDOM> echos out the name of this tmpdir

# usage
if [ "$#" -lt 3 ]; then
    echo "usage: <path/to/probe_util.py> <markersfile> <dir of log2ratio files> [output_dir]"
    exit -1
fi

if [ $4 ]; then
    OUT=$4
else
    # set a random directory if it is not already specified
    OUT=/tmp/cbs_in.$RANDOM
fi

if [ ! -e $OUT ]; then
    # create the directory if necessary
    mkdir $OUT
fi

# full speed ahead
for f in `ls $3/*`;
do
    $1 join_probe_signal $2 $f -o $OUT/`basename $f`
done

echo $OUT
