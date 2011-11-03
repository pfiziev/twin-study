from utils import *
import re
import multiprocessing
from scipy.stats.stats import pearsonr

def read_twdata(fname):
    ids = set()
    sites = {}
    print "reading:", fname
#    for l in open(fname,'r').readlines():
    for l in open(fname,'r'):
        chrNo, nucl, pos, methType, methSubtype, methLevel, methReads, totalReads = filter(None,re.split(r'\s+',l))
        if int(totalReads) >= 4 and methType == 'CG':
            site_id = "%s\t%s" %(chrNo, pos)
            ids.add(site_id)
            sites[site_id] = methLevel
    return sites, ids


def process(tw1, tw2):
    highpriority()
    print "Twin pair: %s %s" % (tw1, tw2)
    tw1_data, tw1_ids = read_twdata(datafiles[tw1])
    elapsed("info %s" % tw1)
    tw2_data, tw2_ids = read_twdata(datafiles[tw2])
    elapsed("info %s" % tw2)

    out = open(os.path.join(DATA_DIR, "%s_%s.sites" % (tw1, tw2)), 'w')
    for site in sorted(tw1_ids & tw2_ids, key = lambda k: map(int, k.split('\t'))):
        out.write("%s\t%s\t%s\n" % (site, tw1_data[site], tw2_data[site]))

    out.close()
    elapsed("%s %s" % (tw1, tw2))


def calc_cc(tw1, (tw1_data, tw1_ids), tw2, (tw2_data, tw2_ids)):
#    highpriority()


    cut = sorted(tw1_ids & tw2_ids, key = lambda k: map(int, k.split('\t')))

    return {'pair' : "random pair: %s %s" % (tw1, tw2),
            'cc' : pearsonr([float(tw1_data[site]) for site in cut], [float(tw2_data[site]) for site in cut]),
            'sites' : len(cut)}

#    print "random pair: %s %s" % (tw1, tw2), pearsonr([float(tw1_data[site]) for site in cut], [float(tw2_data[site]) for site in cut]), 'sites:', len(cut)


if __name__ == '__main__':
#    for tw1, tw2 in twinpairs:
#        multiprocessing.Process(None, process, "%s_%s" % (tw1, tw2), (tw1, tw2)).start()
    data = dict((tw, read_twdata(datafiles[tw])) for tw in twins)


#    multiprocessing.Pool(processes = 4).map(calc_cc, [((tw1, data[tw1]), (tw2, data[tw2])) for tw1, tw2 in random_twin_pairs])
    ccs = [calc_cc(tw1, data[tw1], tw2, data[tw2]) for tw1, tw2 in random_twin_pairs]

    for cc in sorted(ccs, key = lambda x: x['cc'][0]):
        print cc['pair'], 'cc:',cc['cc'][0],'p-value',cc['cc'][1], 'sites:', cc['sites']

    print 'average cc:', sum(cc['cc'][0] for cc in ccs)/len(ccs)


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
