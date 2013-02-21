#!/bin/bash

# runs the function join_probe_signal from probe_utils.py over all the files in a directory
# puts the output in a directory in /tmp/cbs_in.$RANDOM
# echos out the name of this tmpdir

if [ "$#" -ne 3 ]; then
    echo "usage: <path/to/probe_util.py> <markersfile> <dir of log2ratio files>"
    exit -1
fi

OUT=/tmp/cbs_in.$RANDOM
mkdir $OUT
out_dir
for f in `ls $3/*`;
do
    $1 join_probe_signal $2 $f -o $OUT/`basename $f`
done

echo $OUT
