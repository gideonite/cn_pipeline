#!/bin/bash

# takes all the files in the directory (1st argument) and makes a list of them
# in the file given (2nd argument)

for file in `ls $1`
do
    echo $1/$file >> ../CGHPipe/fofs/$2
done
