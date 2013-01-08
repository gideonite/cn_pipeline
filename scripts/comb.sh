#/bin/bash

# Takes an Affymetrix Annotation file (csv), found on the Affymetrix website
# "Support By Product" and turns it into a file for CBS and GISTIC algorithms
# to use

awk '
BEGIN {
    FS=",";
}
{ print $1,"\t",$3,"\t",$4; }
!/---/'
| sed 's/\"//g'


