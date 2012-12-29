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

        hash[probe_name] = [chr, pos]

    sys.stderr.write("DONE!  (skipped %d X/Y probes)\n\n" % skipped_XY)
    return hash

def map_snp_marker_data(hash, snp_lines, diagnostic):
    # takes a hash map like this {'CN_473964': [1, 61808]} and lines from a snp
    # markers file like this
        # header lines
        # AFFX-5Q-123 1267.615
    # and writes to stdout lines like this for input into CBS
    #
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

parser.add_argument('marker_positions', help='file of marker positions, \
        e.g. lines like this SNP_A-1738457    1   328296 ', type=argparse.FileType('r'))

args = parser.parse_args()
args.func(args)

#if __name__ == "__main__":
#    markers_fp = args.marker_positions
#    hash = markers_hash(markers_fp)
#    #print hash
#    markers_fp.close()
#
#    snp_lines = sys.stdin.readlines()
#    #print snp_lines
#
#    if args.unmapped:
#        print 'diagnostic!'
#        print args.map
#
#    if args.map:
#        print 'map away!'
#        #map_snp_marker_data(hash, snp_lines)
#
#    sys.stderr.write("%d unmapped probes\n" % unmapped)
#    sys.stderr.write("%d \n" % len(lines))
