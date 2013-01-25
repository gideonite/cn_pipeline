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
    def __repr__(self):
        return "MarkerSignal(%s, %s)" %(self.name, self.signal)

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

class Level2Data:
    def __init__(self, name, signal):
        self.name = name
        self.signal = signal

def read_cbs_input(cbs_input_f):
    cbs_ins = []

    for line in cbs_input_f.readlines():
        if REF.search(line):
            continue
        line = line.split()
        cbs_in = Level2Data(line[0], line[1])
        cbs_ins.append(cbs_in)
    return cbs_ins

def unmapped_Level2(hash, level2s):
    unmapped = []

    for level2 in level2s:
        try:
            hash[level2.name]
        except KeyError:
            unmapped.append(level2)
    return unmapped

class CbsIn:
    def __init__(self, name, chr, pos, signal):
        self.name = name
        self.chr = chr
        self.pos = pos
        self.signal = signal
    def __repr__(self):
        return "%s\t%s\t%s\t%s" %(self.name, self.chr, self.pos, self.signal)

def join_level2_mp(level2, mp):
    """
    Level2Data, MarkerPos -> CbsIn

    takes a Level2Data object and a MarkerPos object and joins them together,
    forming (and returning!) a CbsIn object
    """
    assert(level2.name == mp.name)
    return CbsIn(mp.name, mp.chr, mp.pos, level2.signal)

def map_level2(hash, level2):
    """
    {probe id : MarkerPos}, Level2Data -> CbsIn
    """
    mp = hash[level2.name]
    return join_level2_mp(level2, mp)
