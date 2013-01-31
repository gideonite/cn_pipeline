#!/bin/bash

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
