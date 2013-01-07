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

# MarkerSignal {{{
class MarkerSignal:
    def __init__(self, name, signal):
        self.name = name
        self.signal = signal
    def __repr__(self):
        return "MarkerSignal(%s, %s)" %(self.name, self.signal)

def read_marker_signal(line):
    #helper function to read in a line
    line = line.strip()
    line = line.split()                          # assume: data has no whitespace
    line = filter(lambda x: x != '', line)

    name, signal = line[0], line[1]
    ms = MarkerSignal(name, signal)
    return ms

def read_marker_signals(_in):
    # reads in a file with lines like this
        # header lines
        # AFFX-5Q-123 1267.615
    # header lines have the string "REF" in them

    out = sys.stdout
    list = []

    sys.stderr.write("reading in marker-signal data...")
    for line in _in.readlines():
        if REF.search(line):
            continue
        marker_signal = read_marker_signal(line)
        list.append(marker_signal)

    sys.stderr.write("DONE!\n")
    return list
#}}}

# {{{ CbsSegment
class CbsSegment:
    def __init__(self, sample, chr, start, stop, seg_mean):
        self.sample = sample
        self.chr = chr
        self.start = start
        self.stop = stop
        self.seg_mean = seg_mean
    def __repr__(self):
        return "CbsSegment(%s, %s, %s, %s, %s)" \
                %(self.sample, self.chr, self.start, self.stop, self.seg_mean)

def read_cbs_segments(_in):
    list = []

    for line in _in.readlines():
        if barcode.search(line):
            continue
        values = line.split()
        cbs_seg = CbsSegment(values[0], values[1], values[2], values[3], values[4])
        list.append(cbs_seg)
    return list
#}}}

class MarkerPosUtil:
    def __init__(self, markerfile):
        self.marker_f = markerfile
        self.chrPos_name_hash = self.make_chrPos_name_hash(markerfile)
        self.locus_hash = None      # initialize

    def make_chrPos_name_hash(self, markerfile):
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

            # not efficient memory-wise but should help safeguard against some
            # errors
            # chr = int(chr)
            # pos = int(pos)

            # 2 way map
            hash[probe_name] = (chr, pos)
            hash[(chr,pos)] = probe_name

        sys.stderr.write("DONE!  (skipped %d X/Y probes)\n\n" % skipped_XY)
        return hash

    def make_locus_hash(self):
        # returns a hash of { chr : position } where position is the locus on the chr

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

    def map_markers(self, markers):
        # map a list of marker names using the hash table
        # and report some stats
        # marker name, e.g. CN_473982

        sys.stderr.write("mapping markers...")
        unmapped = 0
        first_unmapped = ''
        line_no = 0

        for m in markers:
            line_no+=1
            try:
                self.chrPos_name_hash[m]
            except KeyError:
                unmapped += 1
                if first_unmapped == '':
                    first_unmapped = m
                    first_unmapped_line = line_no

        sys.stderr.write("DONE!\n")

        sys.stderr.write("\n%d unmapped probe(s)\n" % unmapped)
        if first_unmapped != '':
            sys.stderr.write("1st unmapped probe: %s (line %d) \n" % (first_unmapped, first_unmapped_line + 1))

    def map_loci(self, loci):
        # map a list of loci using the hash table
        # and report some stats
        # loci, e.g. tuple(1, 82170)

        sys.stderr.write("mapping loci...")
        unmapped = 0
        first_unmapped = ''
        line_no = 0

        for locus in loci:
            try:
                self.chrPos_name_hash[locus]
            except KeyError:
                unmapped += 1
                if first_unmapped == '':
                    first_unmapped = locus
                    first_unmapped_line = line_no
        sys.stderr.write("DONE!\n")

        sys.stderr.write("\n%d unmapped loci\n" % unmapped)
        if first_unmapped != '':
            sys.stderr.write("1st unmapped locus: %s (line %d) \n" % (first_unmapped, first_unmapped_line + 1))

    def nearest_two_probes(self, locus):
        # finds the two probes closest to the locus

        hash = self.get_locus_hash()
        chr = str(locus[0])
        pos = str(locus[1])

        loci = hash[chr]
        loci.sort()     # sorting each time, not sure how expensive this is.

        # these functions are from the pydoc
        # http://docs.python.org/2.7/library/bisect.html#searching-sorted-lists
        def find_le(a, x):
            'Find rightmost value less than or equal to x'
            i = bisect_right(a, x)
            if i:
                return a[i-1]
            raise ValueError
        def find_ge(a, x):
            'Find leftmost item greater than or equal to x'
            i = bisect_left(a, x)
            if i != len(a):
                return a[i]
            raise ValueError

        le = find_le(loci, pos)
        ge = find_ge(loci, pos)

        return (le, ge)

    def distance_to_nearest_probe(self, locus):
        pos = int(locus[1])
        (l, r) = self.nearest_two_probes(locus)
        (l, r) = (int(l), int(r))

        return min(abs(l - pos), abs(r - pos))

if __name__ == "__main__":
    marker_f = open(sys.argv[1])
    util = MarkerPosUtil(marker_f)
    marker_f.close()

    #marker_signals = read_marker_signals(sys.stdin)
    #util.map_markers([ms.name for ms in marker_signals])

    #cbs_segs = read_cbs_segments(sys.stdin)
    #util.map_loci([(seg.chr, seg.start) for seg in cbs_segs])

    #ms.map_marker_signals(read_marker_signals(sys.stdin))
    #ms.map_cbs_segments(read_cbs_segments(sys.stdin))

### {{{ arg parser ###
#parser = argparse.ArgumentParser(description="Utils for dealing with \
#        sequencing probes.  Takes the probes in question via stdin")
#subparsers = parser.add_subparsers()
#
#diagnostic = subparsers.add_parser('diagnostic', help="diagnostic of unmapped probes")
#diagnostic.set_defaults(func=lambda args:
#        map_snp_marker_data(markers_hash(args.marker_positions), \
#            sys.stdin.readlines(), diagnostic=True))
#
#map = subparsers.add_parser('map', help="map away to standard out")
#map.set_defaults(func=lambda args:
#        map_snp_marker_data(markers_hash(args.marker_positions), \
#            sys.stdin.readlines(), diagnostic=False))
#
#map_probes_opt = subparsers.add_parser('map_probes', help="map a list of (newline \
#        deliminited) probes taken from stdin")
#map_probes_opt.set_defaults(func=lambda args:
#        map_probes(markers_hash(args.marker_positions), sys.stdin.readlines()))
#
#map_positions = subparsers.add_parser('map_positions', help="map a list of positions \
#        deliminited) probes taken from stdin")
#map_positions.set_defaults(func=lambda args:
#        map_pos(markers_hash(args.marker_positions), sys.stdin.readlines()))
#
#nearest_probe_opt = subparsers.add_parser('nearest_probe', help="map a list of positions \
#        deliminited) probes taken from stdin")
#nearest_probe_opt.set_defaults(func=lambda args:
#        nearest_probes(chr_pos_hash(args.marker_positions), sys.stdin.readlines()))
#
#parser.add_argument('marker_positions', help='file of marker positions, \
#        e.g. lines like this SNP_A-1738457    1   328296 ', type=argparse.FileType('r'))
#
#args = parser.parse_args()
#args.func(args)
#}}}
