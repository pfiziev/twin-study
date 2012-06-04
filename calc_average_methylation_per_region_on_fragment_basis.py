import gzip
import json
import os
import re
from utils import *

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from numpy.core.multiarray import arange

if __name__ == '__main__':

    anno = { 'All' : set() }
    common_sites = set()
    chroms = set()
    for l in open(os.path.join(ANNO_DIR, 'RRBS_mapable_regions.info.annotated')):
        chr_no, frag_no, frag_start, frag_end, reg_anno = l.split('\t')
        chroms.add(chr_no)
        reg_anno = json.loads(reg_anno)
        for anno_dict in reg_anno:
            if anno_dict['type'] not in anno:
                anno[anno_dict['type']] = set()
            anno[anno_dict['type']].add(frag_no)

        anno['All'].add(frag_no)

    elapsed('annotation')

    fragments = {}
    for tw, fname in datafiles.items():
        fragments[tw] = {}
        for l in open(reg_fname(tw)):

            chr_no, reg_id, start, end, regSites, meth_level = filter(None, re.split(r'\s+', l))
            fragments[tw][reg_id] = float(meth_level)

        elapsed(fname)
    elapsed('reading twin data')

    am =  dict((regType, dict((tw, mean([fragments[tw][frag_id] for frag_id in anno[regType] if frag_id in fragments[tw]]))
                    for tw in twins))
                            for regType in anno)



    elapsed('Average methylation calculation')

    print 'Average Methylation per Region'

    fig = plt.figure()
    ax = plt.subplot(111)
    plt.ylim(0,1)
    plt.xlim(0.5, 10.5)

    x = arange(1,len(twins)+1)
    for rt in sorted(am):
        print rt,'fragments:%d' % len((anno[rt]))  ,'\t', ' '.join(t+':'+str(am[rt][t]) for t in twins)
        if not rt.startswith('repeat') or rt == 'repeat|All':
            ax.plot(x, [am[rt][t] for t in twins], 'o-', label = rt)

    plt.xticks( x,  tuple(twins) )
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    fig.savefig(outd('fragments_methylation_per_region1.png'))

    # plot repeats
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.ylim(0,1)
    plt.xlim(0.5,10.5)


    repTypeNo = 0
    for rt in sorted(am):
        if rt.startswith('repeat'):
            repTypeNo += 1
            if repTypeNo == 8:
                plt.xticks( x ,  tuple(twins) )
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                fig.savefig(outd('fragments_methylation_per_region2.png'))
                fig = plt.figure()
                ax = plt.subplot(111)
                plt.ylim(0,1)
                plt.xlim(0.5,10.5)


            ax.plot(x, [am[rt][t] for t in twins], 'o-', label = rt.split('|')[1])

    plt.xticks( x ,  tuple(twins) )
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    fig.savefig(outd('fragments_methylation_per_region3.png'))

    print '\n\n'
