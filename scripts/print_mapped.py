#!/usr/bin/python

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
    hash = {}
    for mp in mps:
        hash[mp.name] = mp

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

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "usage: <markersfile> <file or dir of level2 data>"
        sys.exit(-1)

    markers_filename = sys.argv[1]
    level2_dir = sys.argv[2]

    # make a hash table out of markerfile
    f = open(markers_filename, 'r')
    mps = read_markerpos(f)
    f.close()
    hash_markerpos = hash_MarkerPos(mps)

    if os.path.isfile(level2_dir):
        level2_filenames = [level2_dir]
    else:
        level2_filenames = os.listdir(level2_dir)

    print 'unmapped probes'
    for filename in level2_filenames:
        basename = os.path.basename(filename)

        f = open(filename, 'r')
        level2s = read_level2_data(f)
        f.close()

        # get all the probes that are unmapped
        unmapped = filter(lambda x: not hash_markerpos.has_key(x.name), level2s)
        print "\t%d\t%s" %(len(unmapped), basename)

        # get all the probes that are mapped,
        # and map them,
        # turning them into cbs input
        mapped = filter(lambda x: hash_markerpos.has_key(x.name), level2s)
        cbs_ins = map(lambda x: map_level2(hash_markerpos, x), mapped)

        # make a cbs_in file
        out = open('cbs_in/' + basename + ".cbs_in", 'w')
        out.truncate()
        out.write(CBS_IN_HEADER + "\n")

        for cbs_in in cbs_ins:
            out_string = repr(cbs_in) + "\n"
            out.write(out_string)

        out.close()
