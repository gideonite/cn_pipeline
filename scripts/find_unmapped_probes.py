#!/usr/bin/python

import sys

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

        probe_name = line[0]
        chr        = line[1].strip()
        pos        = line[2].strip()

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

    if len(sys.argv) == 0:
        print "usage: map_probes <markers_file> stdin"

    # make the hash of markers to chr, pos
    markers_fn = sys.argv[1]
    markers_fp = open(markers_fn)
    hash = markers_hash(markers_fp)
    markers_fp.close()

    lines = sys.stdin.readlines()
    unmapped = 0

    for line in lines:
        line = line.strip()
        line = line.split("\t")

        probe, signal = line[0], line[1]
        try:
            hash[probe][0]
        except KeyError:
            unmapped += 1

    sys.stderr.write("%d unmapped probes\n" % unmapped)
