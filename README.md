Copy Number Pipeline
====================

This is a pipeline for making copy number calls on array data from platforms
like Affymetrix SNP6.0.  Currently, it supports analysis from normalized probe
signals through identification of significantly altered segments.

## Overview of Implementation and Algorithms

Noramlized probe signals are assumed to be in the following format:
`probe_name	signal`

A markersfile performs a join on these normalized probe signals to create a
file of the following format.  It is the input of the algorithm called Circular
Binary Segmentation (CBS):

`probe_name	chromosome	position	signal`

A markersfile looks like this:

`probe_name	chromosome	position`

The output of CBS is a list of segments for every sample.  This list of
segments is then concatenated together and given to the GISTIC algorithm which
does a statistical search for recurrent segments.

## INSTALL

You need to install the algorithms, CBS and GISTIC.

*If you want things to work* out of the box:

1. put your GISTIC installation into `./gistic`

2. put your CBS installation into `./R`

### N.B.

on some systems the MCR might not be able to find the shared library
`libXp.so.6`, for that reason I've included here and modified `LD_LIBRARY_PATH`
per this `lib`, directory.  If while running GISTIC you get an error relating
to this, you may have to deal with it yourself by modifying the env
initialization in `run_gistic.sh`.

### CBS

Look [here](http://www.bioconductor.org/packages/2.11/bioc/html/DNAcopy.html)

### GISTIC

Look [here](http://www.broadinstitute.org/cgi-bin/cancer/publications/pub_paper.cgi?mode=view&paper_id=216&p=t)

### Galaxy

Put the folder containing this README file into the Galaxy `tools` directory.
Append the following to the `toolbox` tag in `tools_conf.xml`:

    <section name="Copy Number Pipeline" id="cn_pipeline">
        <tool file="cn_pipeline/cbs_preprocess.xml" />
        <tool file="cn_pipeline/cbs.xml" />
        <!--<tool file="cn_pipeline/gistic.xml" />-->
    </section>
