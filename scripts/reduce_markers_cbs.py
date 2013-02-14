#!/usr/bin/python

"""
wrapper for the function -- reduce_markers_by_cbs --
this is maninly for preparation for GISTIC
"""

from utils import *
import sys

if len(sys.argv) != 4:
    print "usage: <markersfile> <cbs_file> <output_file>"
    sys.exit(-1)

markers_filename = sys.argv[1]
cbs_file = sys.argv[2]
output_file = sys.argv[3]

# make a hash table out of markerfile
f = open(markers_filename, 'r')
mps = read_markerpos(f)
f.close()
hash_markerpos = hash_MarkerPos(mps)
