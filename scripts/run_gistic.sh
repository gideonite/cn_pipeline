#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "usage: <path/to/gistic/executable> <path/to/MatlabComponentRuntime> <path/to/segfiles> <markersfile> <refgenefile>"
    exit -1
fi

echo --- setting up environment variables ---
mcr_root=$1
export LD_LIBRARY_PATH
export LD_LIBRARY_PATH /home/dresdnerg/galaxy-dist/tools/cn_pipeline/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH $mcr_root/v714/runtime/glnxa64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH $mcr_root/v714/sys/os/glnxa64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH $mcr_root/v714/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH $mcr_root/v714/sys/java/jre/glnxa64/jre/lib/amd64/server:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH $mcr_root/v714/sys/java/jre/glnxa64/jre/lib/amd64:$LD_LIBRARY_PATH
export XAPPLRESDIR $mcr_root/v714/X11/app-defaults

echo --- creating output directory ---
basedir=/tmp/gistic_out.$RANDOM
mkdir $basedir

# merge files
SEG_FILE=$basedir/concated_cbs_out.txt
touch $SEG_FILE
for f in `ls $2/*`;
do
    tail -q -f -n +2 $f >>$SEG_FILE
done

echo --- running GISTIC ---
segfile=$SEG_FILE
markersfile=$3
refgenefile=$4
#alf = $thisdir/examplefiles/arraylistfile.txt
#cnvfile = $thisdir/examplefiles/cnvfile.txt

./gp_gistic2_from_seg -b $basedir -seg $segfile -mk $markersfile -refgene $refgenefile -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90

# input files I'm ignoring for now
#-alf $alf
#-cnv $cnvfile

echo $basedir
