import json
from utils import *
import re
import multiprocessing
from scipy.stats.stats import pearsonr

def read_twdata(fname):
    regions = {}
    print "reading:", fname
    for l in open(fname,'r'):
        chrNo, regId, start, end, sites, methLevel = filter(None, re.split(r'\s+',l))
        if int(sites) >= 3:
            regions[regId] = float(methLevel)
    return regions




def calc_cc(tw1, tw1_data, tw2, tw2_data, filter = None):
    print tw1, tw2

    cut = set(tw1_data) & set(tw2_data)
    if filter is not None:
        cut = cut & filter

    cut = sorted(cut, key = lambda k: map(int, k.split('\t')))


    return {'id' : tw1,
            'pair' : "%s %s" % (tw1, tw2),
            'cc' : pearsonr([float(tw1_data[site]) for site in cut], [float(tw2_data[site]) for site in cut]),
            'sites' : len(cut)}

#    print "random pair: %s %s" % (tw1, tw2), pearsonr([float(tw1_data[site]) for site in cut], [float(tw2_data[site]) for site in cut]), 'sites:', len(cut)


if __name__ == '__main__':

    anno = { 'all' : None }

    for l in open(os.path.join(ANNO_DIR,'RRBS_mapable_regions.info.annotated')):
        _, regId, _, _, regAnno = l.split('\t')
        regAnno = json.loads(regAnno)
        for annoDict in regAnno:
            if annoDict['type'] not in anno:
                anno[annoDict['type']] = set()
            anno[annoDict['type']].add(regId)
    elapsed('annotation')

    data = dict((tw, read_twdata(datafiles[tw]+'.regions')) for tw in twins)
    elapsed('reading twin data')
    ccs = dict((regType, [calc_cc(tw1, data[tw1], tw2, data[tw2], filter = anno[regType]) for tw1, tw2 in twinpairs]) for regType in anno)
    elapsed('calc ccs')


#    ccs = [calc_cc(tw1, data[tw1], tw2, data[tw2]) for tw1, tw2 in twinpairs]
#    ccs = [calc_cc(tw1, data[tw1], tw2, data[tw2]) for tw1, tw2 in random_twin_pairs]
#
    for regType in sorted(ccs):
#        print regType
        sorted_ccs = sorted(ccs[regType], key = lambda x: twin_age[x['id']])
        print "Summary:", regType, "CC:%f\tP-value:%f" % pearsonr([cc['cc'][0] for cc in sorted_ccs],[twin_age[cc['id']] for cc in sorted_ccs])
        for cc in sorted_ccs:
#            print cc['pair'],'\t','age:',twin_age[cc['id']],'\t', 'cc:',cc['cc'][0],'\t','p-value',cc['cc'][1],'\t', 'sites:', cc['sites']
            print regType,'\t', twin_age[cc['id']],'\t', cc['cc'][0],'\t', cc['sites']
        print "\n"
#    print 'average cc:', sum(cc['cc'][0] for cc in ccs)/len(ccs)


#if __name__ == '__main__':
#    for tw1, tw2 in twinpairs:
#        print "Twin pair: %s %s" % (tw1, tw2)
#        tw1_data, tw1_ids = read_twdata(datafiles[tw1])
#        elapsed("info %s" % tw1)
#        tw2_data, tw2_ids = read_twdata(datafiles[tw2])
#        elapsed("info %s" % tw2)
#        out = open(os.path.join(DATA_DIR, "%s_%s.sites" % (tw1, tw2)), 'w')
#        for site in tw1_ids & tw2_ids:
#            out.write("%s, %s\n" % (tw1_data[site], tw2_data[site]))
#
#        out.close()
#        elapsed("done")
