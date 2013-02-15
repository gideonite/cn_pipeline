#!/usr/bin/python

import csv
import sys

columns_aliases = {
        'chrom': 'chr',
        'name': 'probe_id',
        'start': 'pos'
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
    returns a dictionary of key to row.

    key is basically a column header
    """
    return_dict = {}

    for row in rows:
        try:
            val = row[key]
            return_dict[val] = row
        except KeyError:
            print 'no such key:', key
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
        probe_id = s['probe_id']        # let this fail

        try:
            hash[probe_id]['signal'] = s['signal']
        except KeyError:
            unmapped += 1

    sys.stderr.write('unmapped: %d\n' % (unmapped))
    return hash.values()
