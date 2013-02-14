import sys
from utils import *

if len(sys.argv) != 4:
    print "usage: <markersfile> <file or dir of level2 data> <output_dir>"
    sys.exit(-1)

markers_filename = sys.argv[1]
level2_dir = sys.argv[2]
output_dir = sys.argv[3]

# make a hash table out of markerfile
f = open(markers_filename, 'r')
mps = read_markerpos(f)
f.close()
hash_markerpos = hash_MarkerPos(mps)

if os.path.isfile(level2_dir):
    level2_filenames = [level2_dir]
else:
    level2_filenames = os.listdir(level2_dir)

print 'unmapped probes'
for filename in level2_filenames:
    basename = os.path.basename(filename)

    f = open(filename, 'r')
    level2s = read_level2_data(f)
    f.close()

    # get all the probes that are unmapped
    unmapped = filter(lambda x: not hash_markerpos.has_key(x.name), level2s)
    print "\t%d\t%s" %(len(unmapped), basename)

    # get all the probes that are mapped,
    # and map them,
    # turning them into cbs input
    mapped = filter(lambda x: hash_markerpos.has_key(x.name), level2s)
    cbs_ins = map(lambda x: map_level2(hash_markerpos, x), mapped)

    # make a cbs_in file
    out = open(os.path.join(output_dir, basename + ".cbs_in"), 'w')
    out.truncate()
    out.write(CBS_IN_HEADER + "\n")

    for cbs_in in cbs_ins:
        out_string = repr(cbs_in) + "\n"
        out.write(out_string)

    out.close()
