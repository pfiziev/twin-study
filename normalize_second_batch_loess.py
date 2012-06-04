"""
This script normalizes the values of batch 2 according to batch 1.
"""
import json

__author__ = 'pf'

from itertools import izip
from utils import *

if __name__ == '__main__':

    _old_means, new_means_normalized = json.load(open(datad('loess_output_for_fragment_means.json')))

    sites = None

    lines = {}
    data = {}

    for twin_id, fname in datafiles.items():
        fname += '.regions'
        cdata = {}

        for l in open(fname, 'r'):
            buf = l.strip().split()

            frag_id = int(buf[1])
            cdata[frag_id] = float(buf[-1])
            lines[frag_id] = buf

        sites = set(cdata.iterkeys()) if sites is None else sites & set(cdata.iterkeys())
        data[twin_id] = cdata
        elapsed(fname)
        #        print len(sites)

    elapsed('reading data')
    sites = sorted(sites)
    old_twins = sorted([t for t in datafiles if int(t[1:]) <=10])
    new_twins = sorted([t for t in datafiles if int(t[1:]) > 10])
    print old_twins, new_twins

    old_means = {}
    new_means = {}
    for site in sites:
        old_means[site] = mean([data[t][site] for t in old_twins])
        new_means[site] = mean([data[t][site] for t in new_twins])


    elapsed('calculating means')

    for tw in new_twins:

        out = open(datafiles[tw] + '.regions.loess_normalized','w')

        for frag_index, frag_id in enumerate(sites):
            buf = lines[frag_id]

            normalized = data[tw][frag_id] - new_means[frag_id] + new_means_normalized[frag_index]

            # cut values that are below 0 or above 1
            normalized = max(min(normalized, 1), 0)

            buf[-1] = str(normalized)
            out.write('\t'.join(buf)+'\n')

        out.close()
        elapsed('normalized twin: '+ datafiles[tw] + '.regions.loess_normalized')



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

__author__ = 'pf'
