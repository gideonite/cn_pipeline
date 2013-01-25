#!/usr/bin/env python

# Affymetrix SNP6 Level_2 data looks like this,
#
# Hybridization REF AMAZE_p_TCGASNP_b86_87_88_N_GenomeWideSNP_6_A10_735526
# Composite Element REF   Signal
# AFFX-5Q-123 1267.615
# AFFX-5Q-456 314.549
# AFFX-5Q-789 3220.386
# AFFX-5Q-ABC 986.392
# ...

# a probes file looks like this
# SNP_A-1738457    1   328296
# SNP_A-1658232    1   1435232
# ...

import sys
import re
import argparse
from bisect import bisect_left, bisect_right

REF = re.compile("REF")
barcode = re.compile("barcode")

class MarkerSignal: #{{{
    def __init__(self, name, signal):
        self.name = name
        self.signal = signal
    def __repr__(self):
        return "MarkerSignal(%s, %s)" %(self.name, self.signal)
    #}}}

def read_marker_signal(line): #{{{
    #helper function to read in a line
    line = line.strip()
    line = line.split()                          # assume: data has no whitespace
    line = filter(lambda x: x != '', line)

    name, signal = line[0], line[1]
    ms = MarkerSignal(name, signal)
    return ms
#}}}

def read_marker_signals(_in): #{{{
    # reads in a file with lines like this
        # header lines
        # AFFX-5Q-123 1267.615
    # header lines have the string "REF" in them

    list = []
    out = sys.stderr

    out.write("reading in marker-signal data...")
    for line in _in.readlines():
        if REF.search(line):
            continue
        marker_signal = read_marker_signal(line)
        list.append(marker_signal)

    out.write("DONE!\n")
    return list
#}}}

class CbsSegment: #{{{
    def __init__(self, sample, chr, start, stop, seg_mean):
        self.sample = sample
        self.chr = chr
        self.start = start
        self.stop = stop
        self.seg_mean = seg_mean
    def __repr__(self):
        return "CbsSegment(%s, %s, %s, %s, %s)" \
                %(self.sample, self.chr, self.start, self.stop, self.seg_mean)
    def get_start_pos(self):
        return (self.chr, self.start)
    def get_end_pos(self):
        return (self.chr, self.end)
# }}}

def read_cbs_segments(_in): #{{{
    list = []

    out = sys.stderr

    out.write("\nreading cbs segments...")

    for line in _in.readlines():
        if barcode.search(line):
            continue
        values = line.split()
        try:
            cbs_seg = CbsSegment(values[0], values[1], values[2], values[3], values[4])
            list.append(cbs_seg)
        except IndexError:
            raise IndexError(values, "expected sample, chr, start, stop, seg_mean")

    out.write("DONE!\n")

    return list
#}}}

class MarkerPosUtil:
    def __init__(self, markerfile): #{{{
        self.marker_f = markerfile
        self.chrPos_name_hash = self.make_chrPos_name_hash(markerfile)
        self.locus_hash = None      # initialize
        #}}}

    def make_chrPos_name_hash(self, markerfile): #{{{
        # returns a two-way hash (chr, pos) <-> probe_name
        #
        # (someday this may want to be implemented differently, i.e. with a
        # "Marker class")
        sys.stderr.write("hashing probes...")

        lines = markerfile.readlines()

        hash = {}
        skipped_XY = 0

        for line in lines:
            if ("#" in line):
                # skip the first line, is this a dumb way of doing it?
                continue

            line = line.split("\t")

            try:
                probe_name = line[0]
                chr        = line[1].strip()
                pos        = line[2].strip()
            except IndexError:
                print line

            if chr == 'X' or chr == 'Y':
                skipped_XY += 1
                continue

            # 2 way map
            hash[probe_name] = (chr, pos)
            hash[(chr,pos)] = probe_name

        sys.stderr.write("DONE!  (skipped %d X/Y probes)\n" % skipped_XY)
        return hash
    #}}}

    def make_locus_hash(self): #{{{
        # returns a hash of { chr : position } where position is the locus on
        # the chr

        # initialize the map
        hash = {}
        for i in range(1,24):
            if i == 23:
                hash["X/Y"] = []
            else:
                key = str(i)
                hash[key] = []
        # end

        for key in self.chrPos_name_hash.keys():
            if type(key) == tuple:
                hash[key[0]].append(key[1])

        return hash

    def get_locus_hash(self):
        if self.locus_hash == None:
            self.locus_hash = self.make_locus_hash()
        return self.locus_hash
    #}}}

    def unmapped(self, markers): #{{{
        """
        filters through a list of marker_ids OR (chr, pos)
        returning a list of the ones that are not mapped
        (i.e. are not in the dict)
        """

        return filter(lambda m: not self.chrPos_name_hash.has_key(m), markers)
    #}}}

    def nearest_two_probes(self, positions, pos, index): #{{{
        """
        return the two ps in positions that are adjacent to pos.

        index is the index in positions at which to start the search
        """
        # sanity check
        assert(positions[index] < pos)
        positions = sort(positions)

        ps = positions[index:]
        ps_len = len(ps)

        for p_i in xrange(ps_len):
            if pos < ps[p_i]:
                # found a least upper bound
                return (ps[p_i - 1], ps[p_i])

    def distance_to_nearest_probe(self, locus):
        pos = int(locus[1])
        (l, r) = self.nearest_two_probes(locus)
        (l, r) = (int(l), int(r))

        return min(abs(l - pos), abs(r - pos))
    #}}}

def print_unmapped_stats(l, t): #{{{
    """
    counts through the list l and prints out some stats concerning length of l.
    t is a string indicating the type l's elements
    """

    no = len(l)
    sys.stderr.write("\n%d unmapped %s\n" %(no, t))
    if (no < 10):
        for el in l:
            sys.stderr.write("%s\n" %(el,))
#}}}

def unmapped_opt(args): #{{{
    input_type = args.input_file_type
    marker_f = args.marker_file
    input_f = args.input_file

    util = MarkerPosUtil(marker_f)

    if input_type == 'cbs':
        cbs_segs = read_cbs_segments(input_f)
        unmapped = util.unmapped([ seg.get_start_pos() for seg in cbs_segs ])
        print_unmapped_stats(unmapped, "cbs segments")

    elif input_type == 'marker-signal':
        marker_signals = read_marker_signals(input_f)
        unmapped = util.unmapped([ ms.name for ms in marker_signals ])
        print_unmapped_stats(unmapped, "markers")

    input_f.close()
    marker_f.close()
#}}}

def distance_to_nearest_probe_opt(args): #{{{
    marker_f = args.marker_file
    input_f = args.input_file
    out = sys.stderr

    util = MarkerPosUtil(marker_f)

    cbs_segs = read_cbs_segments(input_f)

    seg = cbs_segs[0]
    hash = util.get_locus_hash()
    positions = hash[seg.chr]
    print util.nearest_two_probes(positions, seg.start, 0)

    #print util.nearest_two_probes(range(100), 3, 5)

    #distances = []

    #kout.write("\nfind nearest probes...")
    #kfor seg in cbs_segs:
    #k    5
    #k    #d = util.distance_to_nearest_probe( (seg.chr, seg.start) )

    #k    #if d not in distances:
    #k    #    distances.append(d)
    #k    #    print d
    #kout.write("DONE!\n")
    #}}}

# parser stuff #{{{
parser = argparse.ArgumentParser(description="Utils for dealing with sequencing \
        probes.")
parser.add_argument('--marker_file', '-m', type=argparse.FileType('r'))
parser.add_argument('--input_file', '-i', type=argparse.FileType('r'))
subparsers = parser.add_subparsers()

unmapped = subparsers.add_parser('unmapped', help="print some stats about the \
        number of unmapped probes.")
unmapped.add_argument('--input_file_type', '-t', type=str, help="The type of the \
        input file: {cbs, marker-signal, marker}")
unmapped.set_defaults(func=unmapped_opt)

distance = subparsers.add_parser('distance', help="find the distance to the \
        nearest probe of loci (e.g. chr and positions).  Only takes CBS output as input")
distance.set_defaults(func=distance_to_nearest_probe_opt)

args = parser.parse_args()
args.func(args)
#}}}
