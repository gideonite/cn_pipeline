#!/bin/bash

# converts seg files into the first format for gistic

if [ ! -d $1 ]; then
    echo \"$1\" "is not a directory, sorry"
    exit 1
fi


#tail -q -f -n +2 ls *CBS_out.txt | cut -f1-5,7
