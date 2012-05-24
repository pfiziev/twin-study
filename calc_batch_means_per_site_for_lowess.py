import gzip
from itertools import izip
from utils import *
import json
import cPickle as pickle



__author__ = 'pf'
def hash_site(chr_no, site_pos): return '%s %s' % (chr_no, site_pos)

if __name__ == '__main__':
    sites = set()

    for line in open(os.path.join(ANNO_DIR,'common_sites.annotated')):
        buf = line.split()
        sites.add('%s %s' % (buf[0], buf[1]))

    elapsed('reading common sites')

    data = {}
    for twin_id, fname in datafiles.items():
#        fname += '.regions'
        cdata = {}
#        print fname


        for l in gzip.open(fname, 'r'):

            chrNo, nucl, pos, methType, methSubtype, methLevel, methReads, totalReads = filter(None, re.split(r'\s+',l))
            chrNo = int(chrNo)
            chrNo = 'chr%s' % (str(chrNo) if chrNo <= 22 else ['X','Y','M'][chrNo-23])

            key = '%s %s' % (chrNo, pos)
            if key in sites:
                cdata[key] = float(methLevel)


        data[twin_id] = cdata
        elapsed(fname)
#        print len(sites)

    elapsed('reading data')
    sites = sorted(sites)
    old_twins = [t for t in datafiles if int(t[1:]) <=10]
    new_twins = [t for t in datafiles if int(t[1:]) > 10]
    print old_twins, new_twins

    old_means = []
    new_means = []
    for site in sites:
        old_means.append(mean([data[t][site] for t in old_twins]))
        new_means.append(mean([data[t][site] for t in new_twins]))

    json.dump({'old' : old_means, 'new': new_means}, open(os.path.join(DATA_DIR, 'lowess_means.json'), 'w'))

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

