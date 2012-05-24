import gzip
from itertools import izip
from utils import *
import json
import cPickle as pickle



__author__ = 'pf'
def hash_site(chr_no, site_pos): return '%s %s' % (chr_no, site_pos)

if __name__ == '__main__':
    sites = None


    data = {}
    for twin_id, fname in datafiles.items():
        fname += '.regions'
        cdata = {}
        #        print fname


        for l in open(fname, 'r'):
            _, frag_id, _, _, _, methLevel = l.strip().split()

            cdata[int(frag_id)] = float(methLevel)

        sites = set(cdata.iterkeys()) if sites is None else sites & set(cdata.iterkeys())
        data[twin_id] = cdata
        elapsed(fname)
    #        print len(sites)

    elapsed('reading data')
    sites = sorted(sites)
    old_twins = sorted([t for t in datafiles if int(t[1:]) <=10])
    new_twins = sorted([t for t in datafiles if int(t[1:]) > 10])
    print old_twins, new_twins

    old_means = []
    new_means = []
    for site in sites:
        old_means.append(mean([data[t][site] for t in old_twins]))
        new_means.append(mean([data[t][site] for t in new_twins]))

    json.dump({'old' : old_means, 'new': new_means}, open(os.path.join(DATA_DIR, 'lowess_fragment_means.json'), 'w'))

#    normalized_new_means = lowess(old_means, new_means)
#
#    for twin_id, fname in datafiles.items():
#        if int(twin_id[1:]) > 10:
#            fname += '.regions'
#            out = open(fname+'.normalized', 'w')
#            twin_sites = {}
#
#            for l in open(fname, 'r'):
#                chrNo, regId, start, end, regSites, methLevel = filter(None, re.split(r'\s+',l))
#                regId = int(regId)
#                if regId in sites:
#                    twin_sites[regId] = (l.strip(), float(methLevel))
#
#            for site_index, site in enumerate(sites):
#                line, methLevel = twin_sites[site]
#                out.write('%s\t%f\n' % (line, methLevel - new_means[site_index] + normalized_new_means[site_index]))
#
#            out.close()



## TO PLOT
#d = json.load(open('/mnt/hoffman2/u/home/mcdb/pfiziev/projects/twindata/CGmap/lowess_means.json'))
#
#import numpy as np
#import matplotlib
#import matplotlib.pyplot as plt
#
#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.plot(d['old'], d['new'], '.')
#ax.set_xlabel('old twins')
#ax.set_ylabel('new twins')

