#!/usr/bin/python

import sys
import re
import os

REF = re.compile("REF")

class MarkerPos:
    def __init__(self, name, chr, pos):
        self.name = name
        self.chr = chr
        self.pos = pos
    #def __repr__(self):
    #    return "MarkerSignal(%s, %s)" %(self.name, self.signal)

def read_marker_pos(markerfile_f):
    markerPoses = []

    for line in markerfile_f.readlines():
        line = line.split()

        mp = MarkerPos(line[0], line[1], line[2])
        markerPoses.append(mp)

    return markerPoses

def hashMarkerPos(mps):
    hash = {}
    for mp in mps:
        hash[mp.name] = mp

    return hash

class CBS_in():
    def __init__(self, name, signal):
        self.name = name
        self.signal = signal

def read_cbs_input(cbs_input_f):
    cbs_ins = []

    for line in cbs_input_f.readlines():
        if REF.search(line):
            continue
        line = line.split()
        cbs_in = CBS_in(line[0], line[1])
        cbs_ins.append(cbs_in)
    return cbs_ins

def number_unmapped_cbs_ins(hash, cbs_ins):
    n = 0

    for cbs_in in cbs_ins:
        try:
            hash[cbs_in.name]
        except KeyError:
            n = n + 1
    return n

# hash marker positions
markerfile_name = sys.argv[1]
markerfile_f = open(markerfile_name, 'r')
mps = read_marker_pos(markerfile_f)
markerfile_f.close()

hash = hashMarkerPos(mps)

## read CBS input
#cbsfile_name = sys.argv[2]
#cbsfile_f = open(cbsfile_name, 'r')
#cbs_ins = read_cbs_input(cbsfile_f)
#cbsfile_f.close()

# number unmapped
#print number_unmapped_cbs_ins(hash, cbs_ins)
path = sys.argv[2]
cbs_in_filenames = os.listdir(path)

for cbs_in_filename in cbs_in_filenames:
    print ">>>", cbs_in_filename
    cbsfile_f = open(cbs_in_filename, 'r')
    cbs_ins = read_cbs_input(cbsfile_f)
    cbsfile_f.close()
    print number_unmapped_cbs_ins(hash, cbs_ins)
