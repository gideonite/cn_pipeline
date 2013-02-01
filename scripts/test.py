#!/usr/bin/python

from print_mapped import *

# init some test objects
mps = [MarkerPos('probe 1', '1', '123456'), MarkerPos('probe 2', '2',\
'123456'), MarkerPos('probe 3', '3', '123456'), MarkerPos('probe 4', '4',\
'123456'), MarkerPos('probe 5', '5', '123456')]

hash = hash_MarkerPos(mps)

level2s = [ Level2Data('probe 1', '0.1'), Level2Data('probe 2', '0.2'), \
        Level2Data('probe 3', '0.3'), Level2Data('probe 4', '0.4'),\
        Level2Data('probe 5', '0.5')]

# -- do some tests --

# unmapped_level2
assert([] == unmapped_Level2(hash, level2s))
unmapped = Level2Data('unmapped', '-1')

level2s.append(unmapped)
assert([unmapped] == unmapped_Level2(hash, level2s))

# join_level2_mp
joined = join_level2_mp(Level2Data("probe 1", "0.1"), MarkerPos("probe 1", "23", "12345"))
assert(joined.name == "probe 1")
assert(joined.chr == "23")
assert(joined.pos == "12345")
assert(joined.signal == "0.1")

print "all tests passed"
