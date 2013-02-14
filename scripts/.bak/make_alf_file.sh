#!/bin/bash

if [ ! -f $1 ]; then
    echo \"$1\" "is not a file, sorry"
    exit 1
fi

echo 'array'
cut -f1 $1 | uniq
