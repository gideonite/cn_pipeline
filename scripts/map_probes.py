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

def markers_hash(markerfile):
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

def map_snp_marker_data(hash, snp_lines, diagnostic):
    # takes a hash map like this {'CN_473964': [1, 61808]} and lines from a snp
    # markers file like this
        # header lines
        # AFFX-5Q-123 1267.615
    # header lines are assumed to contain the string "REF" and are ignored
    # based on this criteria
    #
    # diagnostic: boolean. Print things out in format for CBS or just do a diagnostic?

    out = sys.stdout

    REF = re.compile("REF")

    unmapped = 0
    first_unmapped = ''
    line_no = 0
    for line in snp_lines:
        line_no+=1
        if REF.search(line):
            continue
        line = line.strip()
        line = line.split()                          # assume: data has no whitespace
        line = filter(lambda x: x != '', line)

        probe, signal = line[0], line[1]
        try:
            # map the probe
            locus = hash[probe]
            chr = locus[0]
            pos = locus[1]

            if not diagnostic:
                out.write("%s\t%d\t%d\t%s\n" \
                        % (probe, chr, pos, signal))
        except KeyError:
            unmapped += 1
            if first_unmapped == '':
                first_unmapped = probe
                first_unmapped_line = line_no
    sys.stderr.write("\n%d unmapped probe(s)\n" % unmapped)
    if first_unmapped != '':
        sys.stderr.write("1st unmapped probe: %s (line %d) \n" % (first_unmapped, first_unmapped_line))

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

def nearest_probe(markerfile, positions):
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
            positions = hash["X/Y"]
        else:
            positions = hash[chr]

        # no need to duplicate probes
        if positions != [] and pos == positions[-1]:
            continue

        # order is not preserved
        if positions != [] and positions[-1] >= pos:
            sys.stderr.write("Warning: probes are not listed in order of genomic position: ")
            sys.stderr.write("%s %d %d\n" % (chr, positions[-1], pos))
            #sys.exit(1)

        # append the new pos
        positions.append(pos)

    sys.stderr.write("\nDONE hashing!\n")

    for pos in positions:
        print pos


### arg parser ###
parser = argparse.ArgumentParser(description="Utils for dealing with \
        sequencing probes.  Takes the probes in question via stdin")
subparsers = parser.add_subparsers()

diagnostic = subparsers.add_parser('diagnostic', help="diagnostic of unmapped probes")
diagnostic.set_defaults(func=lambda args:
        map_snp_marker_data(markers_hash(args.marker_positions), \
            sys.stdin.readlines(), diagnostic=True))

map = subparsers.add_parser('map', help="map away to standard out")
map.set_defaults(func=lambda args:
        map_snp_marker_data(markers_hash(args.marker_positions), \
            sys.stdin.readlines(), diagnostic=False))

map_probes_opt = subparsers.add_parser('map_probes', help="map a list of (newline \
        deliminited) probes taken from stdin")
map_probes_opt.set_defaults(func=lambda args:
        map_probes(markers_hash(args.marker_positions), sys.stdin.readlines()))

map_positions = subparsers.add_parser('map_positions', help="map a list of positions \
        deliminited) probes taken from stdin")
map_positions.set_defaults(func=lambda args:
        map_pos(markers_hash(args.marker_positions), sys.stdin.readlines()))

nearest_probe_opt = subparsers.add_parser('nearest_probe', help="map a list of positions \
        deliminited) probes taken from stdin")
nearest_probe_opt.set_defaults(func=lambda args:
        nearest_probe(args.marker_positions, sys.stdin.readlines()))

parser.add_argument('marker_positions', help='file of marker positions, \
        e.g. lines like this SNP_A-1738457    1   328296 ', type=argparse.FileType('r'))

args = parser.parse_args()
args.func(args)
