#!/usr/bin/env python

import csv
import sys
import os
import re

columns_aliases = {
        'chrom': 'chr',
        'chromosome': 'chr',
        'name': 'probe_id',
        '#name': 'probe_id',
        'Composite Element REF': 'probe_id',
        'start': 'pos',
        'position': 'pos',
        'Signal': 'signal'
        }

def read_data(filename, delimiter="\t"):
    """
    opens the file, reads it in using DictReader and returns a list of dict
    rows.

    This is basically just a wrapper around DictReader which normalizes column
    names per the columns_aliases
    """
    csv.register_dialect('delimited', delimiter='\t')
    f = open(filename, 'r')

    # ignore a header line that has this : "Hybridization REF".  TODO: If there are
    # more of these then you should make a dict of aliases like above
    first_line = f.readline()
    if not re.search("Hybridization REF", first_line):
        f.seek(0)

    reader = csv.DictReader(f, dialect='delimited')

    rows = []
    for row in reader:
        for key in row:
            try:
                new_key = columns_aliases[key.strip()]
            except KeyError:
                new_key = key
            row[new_key.strip()] = row.pop(key).strip()
        rows.append(row)
    f.close()

    return rows

def make_hash(rows, key):
    """
    returns a dictionary of key to row.  The keys of this new dictionary are
    the values in the column of each row, e.g.
    {'col1':'val1'} -> { 'val1' : {'col1': 'val1'} }
        (if the parameter `key` is 'col1')

    key is basically a column header
    """
    return_dict = {}

    for row in rows:
        try:
            val = row[key]
            return_dict[val] = row
        except KeyError:
            sys.stderr.write('no such key: %s' %(key))
            break

    return return_dict

def join_probe_signal(markers, signals):
    """
    markers : { probe_id, chr, pos }
    signals : { probe_id, signal }

    hashes the markers by probe_id, then joins each signal to it's marker.
    Prints a count of the number of unmapped signals

    returns: list of { probe_id, chr, pos, signal } for input into CBS
    """
    unmapped = 0
    to_return = []

    hash = make_hash(markers, 'probe_id')

    for s in signals:
        try:
            probe_id = s['probe_id']        # let this fail
        except KeyError:
            sys.stderr.write('no column by the name of "probe_id" %s' %(s))
            sys.exit(-1)
        try:
            hash[probe_id]['signal'] = s['signal']
        except KeyError:
            unmapped += 1

    sys.stderr.write('unmapped: %d\n' % (unmapped))

    # filter out the values in the hash map that don't have corresponding
    # signals
    return filter(lambda v: v.has_key('signal'), hash.values())

def print_probe_signal(probe_signals, out):
    """
    prints out a row as tab-delimited, here's an example row:

    {'chr': '1', 'pos': '760188', 'probe_id': 'A_18_P10001394', 'signal': '2.971'}
    """
    header = "%s\t%s\t%s\t%s\n" %( 'probe_id', 'chr', 'pos', 'signal')
    out.write(header)
    for ps in probe_signals:
        try:
            out.write("%s\t%s\t%s\t%s\n" %( ps['probe_id'], ps['chr'], ps['pos'], ps['signal'] ))
        except KeyError:
            sys.stderror('\n%s' %(ps))

            raise KeyError('could not find a column in row', ps, ". Looking for columns ['probe_id', 'chr', 'pos', 'signal']")

def main():
    import argparse

    description = """
    These are some functions for munging probe level data, cbs output, and
    markersfiles
    """
    parser = argparse.ArgumentParser(description=description)

    if len(sys.argv) == 1:
        print parser.parse_args(["-h"])
        sys.exit(0)

    subparsers = parser.add_subparsers()

    join_probe_signal_parser = subparsers.add_parser('join_probe_signal')
    join_probe_signal_parser.add_argument('markersfile', action='store')
    join_probe_signal_parser.add_argument('level_2', action='store')
    join_probe_signal_parser.set_defaults(which="join_probe_signal")
    join_probe_signal_parser.add_argument('-o', '--output_file', action='store')

    filter_markers = subparsers.add_parser('filter_markers')
    filter_markers.add_argument('markersfile', action='store')
    filter_markers.add_argument('cbs_out', action='store')
    filter_markers.set_defaults(which="filter_markers")

    args = parser.parse_args()
    if args.which == 'join_probe_signal':
        level2s = read_data(args.level_2)
        markers = read_data(args.markersfile)
        sys.stderr.write(os.path.basename(args.level_2) + "\t")
        joined = join_probe_signal(markers, level2s)

        if args.output_file:
            out = open(args.output_file, 'w+')
            print_probe_signal(joined, out)
            out.close()
        else:
            print_probe_signal(joined, sys.stdout)

    # for playing around with argparse
    #
    #args = parser.parse_args(["join_probe_signal", "test/cbs_in.txt", "test/agilent_markersfile.txt" ])
    #print args
    #parser.parse_args(["filter_markers", "level2", "markers" ])


if __name__ == "__main__": main()
