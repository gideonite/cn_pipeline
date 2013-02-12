#!/usr/bin/python

import xml.etree.cElementTree as ElementTree
import sys
import re

def parse_xml(xml_file):
    """
    xml_file -> iterator of {'name', 'chr', 'start', 'end'}

    parses out the probe_id, chr, start position, and end position of the probe
    from the xml agilent annotation file
    """

    unmarshallers = {
        'systematic_name': lambda x: dict((
            ('chr', re.sub("chr", "", x.split(":")[0])),
            ('start', x.split(":")[1].split("-")[0]),
            ('end', x.split(":")[1].split("-")[1])
            )),
    }

    loci_test = unmarshallers['systematic_name']("chr7:99480285-99480344")
    assert loci_test['chr'] == '7'
    assert loci_test['start'] == '99480285'
    assert loci_test['end'] == '99480344'

    ignore = ['HsCGHBrightCorner', 'DarkCorner']
    unrecognized = 0

    for event, el in ElementTree.iterparse(xml_file):
        if el.tag == "reporter":
            name= el.attrib['name']
            if name not in ignore:
                try:
                    parsed = unmarshallers['systematic_name'](el.attrib['systematic_name'])
                    parsed['name'] = name
                    yield parsed
                except IndexError as e:
                    #sys.stderr.write("%s not recognized\n" %(name))
                    unrecognized += 1
                el.clear()
    sys.stderr.write("%d unrecognized probes\n" %(unrecognized))

if len(sys.argv) != 2:
    print "usage: <agilent-annotation.xml>. output is printed to stdout"
    sys.exit(-1)

xml_filename = sys.argv[1]

xml_f = open(xml_filename, 'r')

print 'name', '\t' 'chr', '\t', 'start', '\t', 'end', '\t'
for x in parse_xml(xml_f):
    print x['name'], '\t', x['chr'], '\t', x['start'], '\t', x['end']

xml_f.close()
