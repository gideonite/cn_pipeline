#!/usr/bin/python

from probe_utils import *
import sys

def read_data_test():
    print "--read_data test--"
    cbs_outs = read_data('test/agilent_cbs.txt')
    markers = read_data('test/agilent_markersfile.txt')
    log2ratios = read_data('test/log2ratio.txt')
    tcga_log2ratios = read_data('test/tcga-log2ratio.txt')

    # test for field/col names ONLY, no testing for values
    # TODO: Do you want to test for values?  Maybe? =)
    print '\ttesting for field/col names (ONLY! no testing done on values)'
    cbs_outs[0]['seg.mean']
    cbs_outs[0]['loc.start']
    cbs_outs[0]['loc.end']
    cbs_outs[0]['chr']
    cbs_outs[0]['num.mark']
    cbs_outs[0]['ID']

    markers[0]['chr']
    markers[0]['pos']
    markers[0]['probe_id']

    log2ratios[0]['signal']
    log2ratios[0]['probe_id']

    tcga_log2ratios[0]['signal']
    tcga_log2ratios[0]['probe_id']

def make_hash_test():
    print "\n--make_hash test--"

    row = dict(( ('col1', 'val1'), ('col2','val2') ))
    assert(make_hash([row], 'col1')['val1'] == row)

def print_probe_signal_test():
    print "\n--print_probe_signal_test--"
    try:
        print_probe_signal([{'pos': '760188', 'probe_id': 'A_18_P10001394', 'signal': '2.971'}], sys.stdout)
    except KeyError, e:
        print "\t--print_probe_signal passed--"
        #print e

def join_probe_signal_test():
    print "\n--join_probe_signal_test--"
    markers = [ dict(( ('probe_id', 'ABC_123'), ('chr', 'X'), ('pos', '1') )),
            dict(( ('probe_id', 'ABC_1111'), ('chr', 'X'), ('pos', '1') )) ]
    signals = [ dict(( ('probe_id', 'ABC_123'), ('signal', 'whahoooooooBOOM!') )),
                dict(( ('probe_id', 'ABC_345'), ('signal', 'whahoooooooBOOM!') )) ]

    joined = join_probe_signal(markers, signals)

    assert len(joined) == 1

    assert joined[0]['probe_id'] == 'ABC_123'
    assert joined[0]['chr'] == 'X'
    assert joined[0]['signal'] == 'whahoooooooBOOM!'

def sort_cbs_out_test():
    print "\n--sort_cbs_out_test--"

    zero   = dict(( ('ID', 'id-1'), ('chr', 1), ('loc.start', 0), ('loc.end', 25611452), ('num.mark', 3), ('seg.mean', 5.4943)))
    one    = dict(( ('ID', 'id-2'), ('chr', 1), ('loc.start', 1), ('loc.end', 25611452), ('num.mark', 3), ('seg.mean', 5.4943)))
    two    = dict(( ('ID', 'id-3'), ('chr', 2), ('loc.start', 0), ('loc.end', 25611452), ('num.mark', 3), ('seg.mean', 5.4943)))
    three  = dict(( ('ID', 'id-4'), ('chr', 2), ('loc.start', 1), ('loc.end', 25611452), ('num.mark', 3), ('seg.mean', 5.4943)))
    four   = dict(( ('ID', 'id-5'), ('chr', 'X'), ('loc.start', 0), ('loc.end', 25611452), ('num.mark', 3), ('seg.mean', 5.4943)))
    five   = dict(( ('ID', 'id-6'), ('chr', 'Y'), ('loc.start', 1), ('loc.end', 25611452), ('num.mark', 3), ('seg.mean', 5.4943)))

    cbs_outs = [ two, three, one, zero, five, four ]

    sort_cbs_out(cbs_outs)

    assert cbs_outs[0] == zero
    assert cbs_outs[1] == one
    assert cbs_outs[2] == two
    assert cbs_outs[3] == three
    assert cbs_outs[4] == four
    assert cbs_outs[5] == five

def format_cbs_out_test():
    print "\n--format_cbs_out_test--"

    ex = dict(( ('ID', 'id-2'), ('chr', 1), ('loc.start', 1), ('loc.end', 25611452), ('num.mark', 3), ('seg.mean', 5.4943) ))

    print format_cbs_out('')
    assert format_cbs_out('') == "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean"
    assert format_cbs_out(ex) == "id-2\t1\t1\t25611452\t3\t5.4943"

read_data_test()
make_hash_test()
print_probe_signal_test()
join_probe_signal_test()
sort_cbs_out_test()
format_cbs_out_test()

print "\nall tests passed"
