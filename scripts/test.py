#!/usr/bin/python

from utils import *

print "read_data test"
cbs_outs = read_data('test/agilent_cbs.txt')
markers = read_data('test/agilent_markersfile.txt')
level_2s = read_data('test/cbs_in.txt')

#print cbs_outs
#print markers
#print level_2s
#
#print make_hash(cbs_outs, 'name')
#print make_hash(markers, 'name')
#print make_hash(level_2s, 'name')

joined = join_probe_signal(markers, level_2s)
#for x in joined:
#    print x
print_probe_signal(joined)

print "\nall tests passed"
