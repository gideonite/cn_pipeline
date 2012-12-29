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
import argparse

def markers_hash(markerfile):
    sys.stderr.write("hashing probes...\n")

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

    sys.stderr.write("...DONE!  (skipped %d X/Y probes)\n" % skipped_XY)
    return hash

if __name__ == "__main__":

    # make an arg parser
    parser = argparse.ArgumentParser(description="Utils for dealing with \
            sequencing probes.  Takes the probes in question via stdin")

    # some options
    parser.add_argument('unmapped', help="diagnostic of unmapped probes")
    parser.add_argument('map', help="map away")
    parser.add_argument('--to', help="file to map probes into", type=argparse.FileType('r'))
    parser.add_argument('markers_file', help='file of probes', type=argparse.FileType('r'))

    args = parser.parse_args()

    ## make the hash of markers to chr, pos
    #markers_fn = sys.argv[1]
    #markers_fp = open(markers_fn)
    #hash = markers_hash(markers_fp)
    #markers_fp.close()

    #lines = sys.stdin.readlines()
    #unmapped = 0

    #for line in lines:
    #    line = line.strip()
    #    line = line.split("\t")

    #    probe, signal = line[0], line[1]
    #    try:
    #        #sys.stdout.write("%s\t%d\t%d\t%s\n" % (probe, hash[probe][0], hash[probe][1], signal))     # print map string
    #        hash[probe][0]                      # map the probe
    #    except KeyError:
    #        unmapped += 1

    #sys.stderr.write("%d unmapped probes\n" % unmapped)
    #sys.stderr.write("%d \n" % len(lines))
