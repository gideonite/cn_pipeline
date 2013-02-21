#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "usage: <path/to/run_cbs.r> <dir of cbs_in files> <log_file>"
    exit -1
fi

OUT=/tmp/cbs_out.$RANDOM
mkdir $OUT

for f in `ls $2/*`;
do
    R CMD BATCH --no-save "--args $f $OUT/`basename $f`" $1 $3
    #`basename $f`.log
done

echo $OUT
