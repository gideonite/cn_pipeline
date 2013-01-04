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

REF = re.compile("REF")

class MarkerSignal:
    def __init__(self, name, signal):
        self.name = name
        self.signal = signal
    def __repr__(self):
        return "MarkerSignal(%s, %s)" %(self.name, self.signal)

#{{{
#def map_snp_marker_data(hash, snp_lines, diagnostic):
#    # takes a hash map like this {'CN_473964': [1, 61808]} and lines from a snp
#    # markers file like this
#        # header lines
#        # AFFX-5Q-123 1267.615
#    # header lines are assumed to contain the string "REF" and are ignored
#    # based on this criteria
#    #
#    # diagnostic: boolean. Print things out in format for CBS or just do a diagnostic?
#
#    out = sys.stdout
#
#    REF = re.compile("REF")
#
#    unmapped = 0
#    first_unmapped = ''
#    line_no = 0
#    for line in snp_lines:
#        line_no+=1
#        if REF.search(line):
#            continue
#        line = line.strip()
#        line = line.split()                          # assume: data has no whitespace
#        line = filter(lambda x: x != '', line)
#
#        probe, signal = line[0], line[1]
#        try:
#            # map the probe
#            locus = hash[probe]
#            chr = locus[0]
#            pos = locus[1]
#
#            if not diagnostic:
#                out.write("%s\t%d\t%d\t%s\n" \
#                        % (probe, chr, pos, signal))
#        except KeyError:
#            unmapped += 1
#            if first_unmapped == '':
#                first_unmapped = probe
#                first_unmapped_line = line_no
#    sys.stderr.write("\n%d unmapped probe(s)\n" % unmapped)
#    if first_unmapped != '':
#        sys.stderr.write("1st unmapped probe: %s (line %d) \n" % (first_unmapped, first_unmapped_line))
#}}}

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

    unmapped = 0
    first_unmapped = ''
    line_no = 0
    list = []

    for line in _in.readlines():
        line_no+=1
        if REF.search(line):
            continue
        markerSignal = read_marker_signal(line)
        list.append(markerSignal)

    return list

class CbsSegment:
    def __init__(self, sample, chr, start, stop, num_mark, seg_mean):
        self.sample = sample
        self.chr = chr
        self.start = start
        self.stop = stop
        self.num_mark = num_mark
        self.seg_mean = seg_mean

class MarkerPosUtil:
    def __init__(self, markerfile):
        self.marker_f = markerfile
        self.hash = self.make_hash(markerfile)

    def make_hash(self, markerfile):
        # someday this may want to be implemented differently, i.e. with a
        # "Marker class"
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
            chr = int(chr)
            pos = int(pos)

            # 2 way map
            hash[probe_name] = (chr, pos)
            hash[(chr,pos)] = probe_name

        sys.stderr.write("DONE!  (skipped %d X/Y probes)\n\n" % skipped_XY)
        return hash

def map_probes(hash, probes):
    # takes a list of probes and looks for unmapped probes

    out = sys.stderr
    first_unmapped_probe, first_unmapped_probe_line = '', -1
    unmapped = 0

    i = 0
    for probe in probes:
        i += 1
        try:
            hash[probe]
        except KeyError:
            if first_unmapped_probe == '':
                first_unmapped_probe = probe
                first_unmapped_line = i
            unmapped += 1
    out.write("\n%d unmapped probe(s)\n" % unmapped)

    if first_unmapped_probe != '':
        out.write("1st unmapped probe: %s (line %d) \n" \
                % (first_unmapped_probe, first_unmapped_line))

def map_pos(hash, positions):
    # takes a list of positions and looks for unmapped positions
    # positions are a list of strings like this
        # chr   position
        # 1     1234
    out = sys.stderr

    unmapped = 0
    for position in positions:
        pos = tuple(position.split())
        try:
            hash[pos]
        except KeyError:
            unmapped += 1
    out.write("\n%d unmapped position(s)\n" % unmapped)

def chr_pos_hash(markerfile):
    sys.stderr.write("hashing probes...")

    lines = markerfile.readlines()

    # initialize hash map of { chr : position }
    hash = {}
    for i in range(1,24):
        if i == 23:
            hash["X/Y"] = []
        else:
            key = str(i)
            hash[key] = []

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

        pos = int(pos)

        # get the proper mapping
        if chr == 'X' or chr == 'Y':
            loci = hash["X/Y"]
        else:
            loci = hash[chr]

        # no need to duplicate probes
        if loci != [] and pos == loci[-1]:
            continue

        # order is not preserved
        if loci != [] and loci[-1] >= pos:
            sys.stderr.write("Warning: probes are not listed in order of genomic position: ")
            sys.stderr.write("%s %d %d\n" % (chr, loci[-1], pos))
            #sys.exit(1)

        # append the new pos
        loci.append(pos)

    sys.stderr.write("\nDONE hashing!\n")
    return hash

def nearest_probe(hash, lookup):
    # lookup may be a position or a probe name, anything that might be found in
    # the hash
    print lookup
    print hash[lookup]

def nearest_probes(hash, lookups):

    REF = re.compile("REF")

    print nearest_probe(hash, lookups[3])


if __name__ == "__main__":
    marker_f = open(sys.argv[1])
    ms = MarkerPosUtil(marker_f)
    marker_f.close()

    print read_marker_signals(sys.stdin)


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
