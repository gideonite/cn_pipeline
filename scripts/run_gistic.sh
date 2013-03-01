#!/bin/bash

if [ "$#" -ne 5 ]; then
    echo "usage: <path/to/gistic/executable> <path/to/MatlabComponentRuntime> <path/to/segfiles> <markersfile> <refgenefile>"
    exit -1
fi

EXEC=$1
MCR=$2
SEG_FILES=$3
MARKERSFILE=$4
REFGENEFILE=$5

echo --- setting up environment variables ---
mcr_root=$MCR
export LD_LIBRARY_PATH=/home/dresdnerg/galaxy-dist/tools/cn_pipeline/lib
export LD_LIBRARY_PATH=$mcr_root/v714/runtime/glnxa64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$mcr_root/v714/sys/os/glnxa64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$mcr_root/v714/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$mcr_root/v714/sys/java/jre/glnxa64/jre/lib/amd64/server:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$mcr_root/v714/sys/java/jre/glnxa64/jre/lib/amd64:$LD_LIBRARY_PATH
export XAPPLRESDIR=$mcr_root/v714/X11/app-defaults

echo $LD_LIBRARY_PATH

echo --- creating output directory ---
basedir=/ifs/e63data/sander-lab/dresdnerg/galaxy-tmp/gistic_out.$RANDOM
mkdir $basedir

# merge files
echo --- merging files ---
MERGED_SEGS=/ifs/e63data/sander-lab/dresdnerg/galaxy-tmp/concated_cbs_out.$RANDOM.txt

CURR_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
for f in `ls $SEG_FILES/*`;
do
    #quick and dirty
    #to sort: replace,  X -> 23, Y ->24, then sort, then replace them back
    #tail -q -n +2 $f \
    #| awk '{ if ($2=="X") {$2=="23"} }1' | awk '{ if ($2=="Y") {$2=="24"} }1' \
    #| sort -k 2,2n -k 3,3n \
    #| awk '{ if ($2=="23") {$2=="X"} }1' | awk '{ if ($2=="24") {$2=="Y"} }1' \
    #./
    # >>$MERGED_SEGS

    $CURR_DIR/probe_utils.py sort_cbs $f | tail -q -n +2 >>$MERGED_SEGS
done

echo $MERGED_SEGS

echo --- running GISTIC ---
#alf = $thisdir/examplefiles/arraylistfile.txt
#cnvfile = $thisdir/examplefiles/cnvfile.txt

$EXEC -b $basedir -seg $MERGED_SEGS -mk $MARKERSFILE -refgene $REFGENEFILE -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90

mv $MERGED_SEGS $basedir

# input files I'm ignoring for now
#-alf $alf
#-cnv $cnvfile

echo $basedir
