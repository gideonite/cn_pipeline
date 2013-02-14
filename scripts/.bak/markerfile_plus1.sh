#!/bin/bash

# adds 1 to probe locations

if [ ! -f $1 ]; then
    echo \"$1\" "is not a file, sorry"
    exit 1
fi

if [ -f $2 ]; then
    echo \"$2\" "exists. aborting"
fi

echo "##########THE FOLLOWING PROBES ARE SKIPPED#######" >&2
awk '!/random/ { print $5, substr($2, 4), $3 + 1 }
    /random/ { print $0 >"/dev/stderr"; }' $1
    #$1 >$2

echo "##########THE ABOVE PROBES WERE SKIPPED#######" >&2
