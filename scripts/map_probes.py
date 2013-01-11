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

class MarkerPosUtil: #{{{
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

            # 2 way map
            hash[probe_name] = (chr, pos)
            hash[(chr,pos)] = probe_name

        sys.stderr.write("DONE!  (skipped %d X/Y probes)\n" % skipped_XY)
        return hash

    def make_locus_hash(self):
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

    def unmapped_markers(self, markers):
        """
        filters through a list of marker ids,
        returning a list of the ones that are not mapped
        (i.e. do not have a corresponding chromosome and position in the hash
        table)
        """

        return filter(lambda m: not self.chrPos_name_hash.has_key(m), markers)

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

    def nearest_two_probes(self, locus, index):
        # finds the two probes closest to the locus

        # loops over the locus' in the chromosome, popping as it goes

        hash = self.get_locus_hash()
        chr = str(locus[0])
        pos = str(locus[1])

        loci = hash[chr]
        loci.sort()     # sorting each time, not sure how expensive this is.

        #for locus in loci

        # these functions are from the pydoc
        # http://docs.python.org/2.7/library/bisect.html#searching-sorted-lists
        def find_le(a, x):
            'Find rightmost value less than or equal to x'
            i = bisect_right(a, x)
            if i:
                return a[i-1]
            raise ValueError("chromsome probably not found")
        def find_ge(a, x):
            'Find leftmost item greater than or equal to x'
            i = bisect_left(a, x)
            if i != len(a):
                return a[i]
            raise ValueError("chromsome probably not found")

        le = find_le(loci, pos)
        ge = find_ge(loci, pos)

        return (le, ge)

    def distance_to_nearest_probe(self, locus):
        pos = int(locus[1])
        (l, r) = self.nearest_two_probes(locus)
        (l, r) = (int(l), int(r))

        return min(abs(l - pos), abs(r - pos))
#}}}

def print_unmapped(l, t):
    """
    counts through the list l and prints out some stats concerning length of l.
    t is a string indicating the type l's elements
    """

    no = len(l)
    sys.stderr.write("\n%d unmapped %s\n" %(no, t))
    if (no < 10):
        for el in l:
            sys.stderr.write("%s\n" %(el))

def unmapped_opt(args):
    input_type = args.input_file_type
    marker_f = args.marker_file
    input_f = args.input_file

    util = MarkerPosUtil(marker_f)

    if input_type == 'cbs':
        cbs_segs = read_cbs_segments(input_f)
        util.map_loci([ (seg.chr, seg.start) for seg in cbs_segs])
    elif input_type == 'marker-signal':
        marker_signals = read_marker_signals(input_f)
        unmapped = util.unmapped_markers([ ms.name for ms in marker_signals ])
        print_unmapped(unmapped, "markers")

    input_f.close()
    marker_f.close()

def distance_to_nearest_probe_opt(args):
    marker_f = args.marker_file
    input_f = args.input_file
    out = sys.stderr

    util = MarkerPosUtil(marker_f)

    cbs_segs = read_cbs_segments(input_f)

    distances = []

    out.write("\nfind nearest probes...")
    for seg in cbs_segs:
        5
        #d = util.distance_to_nearest_probe( (seg.chr, seg.start) )

        #if d not in distances:
        #    distances.append(d)
        #    print d
    out.write("DONE!\n")

# parser stuff
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
