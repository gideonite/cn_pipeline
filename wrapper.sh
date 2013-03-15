#!/bin/bash

# basically experimenting with SGE_TASK_ID and keep it for the sake of the
# record, perhaps at some point this knowledge will be useful for writing other
# scripts to submit to cluster

i=0
#SGE_TASK_ID=3      # test

for f in `ls | sort`;
do
    let i=i+1;

    if [ $i == $SGE_TASK_ID ];
    then
        echo $f
    fi
done
