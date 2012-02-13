import json
import os
import re
from utils import *
import matplotlib.pyplot as plt
from numpy.core.multiarray import arange

if __name__ == '__main__':

    anno = { 'All' : set() }
    common_sites = set()
    chroms = set()
    for l in open(os.path.join(ANNO_DIR,'common_sites.annotated')):
        chr_no, site_pos, regAnno = l.split('\t')
        chroms.add(chr_no)
        regAnno = json.loads(regAnno)
        for annoDict in regAnno:
            if annoDict['type'] not in anno:
                anno[annoDict['type']] = set()
            anno[annoDict['type']].add((chr_no, site_pos))

        anno['All'].add((chr_no, site_pos))
        
    elapsed('annotation')

    sites = {}
    for fname in os.listdir(DATA_DIR):
        if fname.endswith('.CGmap'):
            tw = fname.split('_')[0]
            sites[tw] = dict((chrom, {}) for chrom in chroms)
            for l in open(os.path.join(DATA_DIR, fname), 'r'):
                chrNo, nucl, pos, methType, methSubtype, methLevel, methReads, totalReads = filter(None,re.split(r'\s+',l))
                chrNo = int(chrNo)
                if int(totalReads) >= 4 and methType == 'CG':
                    sites[tw]['chr%s' % (str(chrNo) if chrNo <= 22 else ['X','Y','M'][chrNo-23])][pos] = float(methLevel)
            elapsed(fname)
    elapsed('reading twin data')

    am =  dict((regType, dict((tw, sum(sites[tw][chrNo][pos] for chrNo, pos in anno[regType])/len(anno[regType]))
                                        for tw in twins))
                                            for regType in anno)



    elapsed('Average methylation calculation')

    print 'Average Methylation per Region'

    fig = plt.figure()
    ax = plt.subplot(111)
    plt.ylim(0,1)
    plt.xlim(0.5, 10.5)

    x = arange(1,11)
    for rt in sorted(am):
        print rt,'sites:%d' % len((anno[rt]))  ,'\t', ' '.join(t+':'+str(am[rt][t]) for t in twins)
        if not rt.startswith('repeat') or rt == 'repeat|All':
            ax.plot(x, [am[rt][t] for t in twins], 'o-', label = rt)

    plt.xticks( x,  tuple(twins) )
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    fig.savefig('sites_methylation_per_region1.png')

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
                fig.savefig('sites_methylation_per_region2.png')
                fig = plt.figure()
                ax = plt.subplot(111)
                plt.ylim(0,1)
                plt.xlim(0.5,10.5)


            ax.plot(x, [am[rt][t] for t in twins], 'o-', label = rt.split('|')[1])

    plt.xticks( x ,  tuple(twins) )
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    fig.savefig('sites_methylation_per_region3.png')

    print '\n\n'
