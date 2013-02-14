#!/usr/bin/python

"""
Utils for dealing with markerfiles, cbs input, and cbs output.

Contains:
- classes for each data type
- functions for reading in data from files
- functions for doing various filtering of data

Gideon Dresdner     <dresdnerg@cbio.mskcc.org>
"""

import sys
import re
import os

REF = re.compile("REF")
CBS_IN_HEADER="probe_id\tchr\tpos\tsignal"      # important : this is synched with run_cbs.r

class MarkerPos:
    def __init__(self, name, chr, pos):
        self.name = name
        self.chr = chr
        self.pos = pos
    def __repr__(self):
        return "MarkerPos(%s, %s, %s)" %(self.name, self.chr, self.pos)

def read_markerpos(markerfile_f):
    markerPoses = []

    for line in markerfile_f.readlines():
        line = line.split()

        mp = MarkerPos(line[0], line[1], line[2])
        markerPoses.append(mp)

    return markerPoses

def hash_MarkerPos(mps):
    """
    returns {name, (chr, pos) -> MarkerPos}
    a hashmap with name keys and [chr, pos] keys that both map to MarkerPos
    """
    hash = {}
    for mp in mps:
        hash[mp.name] = mp
        hash[ (chr, mp.pos) ] = mp

    return hash

class Level2Data:
    def __init__(self, name, signal):
        self.name = name
        self.signal = signal
    def __repr__(self):
        return "Level2Data(%s, %s)" %(self.name, self.signal)
    def __eq__(self, other):
        return self.name == other.name and self.signal == other.signal

def read_level2_data(cbs_input_f):
    """
    file -> [Level2Data]
    """

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
    """
    CBS input format
    """
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

class CbsOut:
    """
    CBS output format
    """
    def __init__(self, sample_id, chr, seg_start, seg_end, num_markers, seg_mean):
        self.sample_id = sample_id
        self.chr = chr
        self.seg_start = seg_start
        self.seg_end = seg_end
        self.num_markers = num_markers
        self.seg_mean = seg_mean
    def __repr__(self):
        return "%s\t%s\t%s\t%s\t%s\t%s" %(self.sample_id, self.chr,\
                self.seg_start, self.seg_end, self.num_markers, self.seg_mean)

def read_cbs_out(filename):
    """
    reads in a cbs file
    """
    cbs_outs = []

    ID = re.compile("ID")
    for line in filename.readlines():
        if not ID.search(line):
            line = line.split()
            cbs_out = CbsOut(line[0], line[1], line[2], line[3], line[4], line[5])
            cbs_outs.append(cbs_out)

    return cbs_outs

def reduce_markers_by_cbs(hash, cbs_outs):
    """
    {[chr, pos] -> MarkerPos} , [CbsIn] --> [MarkerPos]
    returns all the MarkerPos that have a corresponding CBS segment
    """
    mps = []
    for cbs_out in cbs_outs:
        try:
            mp = hash[ (cbs_out.chr, cbs_out.seg_start) ]
            mps.append(mp)
        except KeyError:
            print "unmapped segment", cbs_out
            sys.exit(-1)
    return mps

