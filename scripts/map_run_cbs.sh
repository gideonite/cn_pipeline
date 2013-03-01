#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "usage: <path/to/run_cbs.r> <path/to/DNAcopy> <dir of cbs_in files> <log_file>"
    exit -1
fi

OUT=/ifs/e63data/sander-lab/dresdnerg/galaxy-tmp/cbs_out.$RANDOM
mkdir $OUT

R=/srv/opt/R/bin/R
for f in `ls $3/*`;
do
    #R CMD BATCH --no-save "--args path/to/DNAcopy in_file out_file" /path/to/run_cbs.r cbs.log
    $R CMD BATCH --no-save "--args $2 $f $OUT/`basename $f`" $1 $4
    #`basename $f`.log
done

echo $OUT
