Copy Number Pipeline
--------------------

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
segments is then concatenated together and given to the GISTIC algorithm for a
statistical search of recurrent segments.

## INSTALL

### Galaxy

Put the folder containing this README file into the Galaxy `tools` directory.
Append the following to the `toolbox` tag in `tools_conf.xml`:

    <section name="Copy Number Pipeline" id="cn_pipeline">
        <tool file="cn_pipeline/cbs_preprocess.xml" />
        <tool file="cn_pipeline/cbs.xml" />
        <!--<tool file="cn_pipeline/gistic.xml" />-->
    </section>
