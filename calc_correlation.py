import matplotlib
matplotlib.use('Agg')

import json
from numpy.core.multiarray import arange
import pylab
from utils import *
import re
import multiprocessing


from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt

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

    anno = { 'All' : None }

    for l in open(os.path.join(ANNO_DIR,'RRBS_mapable_regions.info.annotated')):
        _, regId, _, _, regAnno = l.split('\t')
        regAnno = json.loads(regAnno)
        for annoDict in regAnno:
            if annoDict['type'] not in anno:
                anno[annoDict['type']] = set()
            anno[annoDict['type']].add(regId)
    elapsed('annotation')

    data = dict((tw, read_twdata(reg_fname(tw))) for tw in twins)
    elapsed('reading twin data')


    # plot average methylation per region

    cut = reduce(lambda x,y: x&y, (set(data[t]) for t in data ))
    def calc_am(sites, cut, filter = None):
        cut = (cut & filter if filter is not None else cut)
        return sum(sites[f] for f in cut)/len(cut)

    am = dict((regType, dict((tw, calc_am(data[tw], cut, filter = anno[regType])) for tw in twins)) for regType in anno)



    print 'Average Methylation per Region'

    fig = plt.figure()
    ax = plt.subplot(111)
    plt.ylim(0,1)
    plt.xlim(0.5, 10.5)

    x = arange(1,len(twins) + 1)
    for rt in sorted(am):
        print rt,'fragments:%d' % len((cut & anno[rt]) if anno[rt] is not None else cut)  ,'\t', ' '.join(str(am[rt][t]) for t in twins)
        if not rt.startswith('repeat') or rt == 'repeat|All':
            ax.plot(x, [am[rt][t] for t in twins], 'o-', label = rt)

    plt.xticks( x,  tuple(twins) )
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    fig.savefig(outd('methylation_per_region1.png'))

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
                fig.savefig(outd('methylation_per_region2.png'))
                fig = plt.figure()
                ax = plt.subplot(111)
                plt.ylim(0,1)
                plt.xlim(0.5,10.5)


            ax.plot(x, [am[rt][t] for t in twins], 'o-', label = rt.split('|')[1])

    plt.xticks( x ,  tuple(twins) )
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    fig.savefig(outd('methylation_per_region3.png'))

    print '\n\n'


    # calculate and plot correlation

    ccs = dict((regType, [calc_cc(tw1, data[tw1], tw2, data[tw2], filter = anno[regType]) for tw1, tw2 in twinpairs]) for regType in anno)
    elapsed('calc ccs')

    summary_cc = {}
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.xlim(xmin = 20,xmax = 56)

    for regType in sorted(ccs):
#        print regType
        sorted_ccs = sorted(ccs[regType], key = lambda x: twin_age[x['id']])
        summary_cc[regType] = pearsonr([cc['cc'][0] for cc in sorted_ccs],[twin_age[cc['id']] for cc in sorted_ccs])
        print "Summary:", regType, "CC:%f\tP-value:%f" % summary_cc[regType]
        for cc in sorted_ccs:
#            print cc['pair'],'\t','age:',twin_age[cc['id']],'\t', 'cc:',cc['cc'][0],'\t','p-value',cc['cc'][1],'\t', 'sites:', cc['sites']
            print regType,'\t', twin_age[cc['id']],'\t', cc['cc'][0],'\t', cc['sites']

        if not regType.startswith('repeat') or regType == 'repeat|All':
            ax.plot([twin_age[cc['id']] for cc in sorted_ccs],[cc['cc'][0] for cc in sorted_ccs], 'o-', label = regType)


        print "\n"
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig.savefig(outd('correlation_per_region_type.png'))

    fig = plt.figure()
    ax = plt.subplot(111)

    plt.xlim(xmin = 20, xmax = 56)

    repTypeNo = 0
    for regType in sorted(ccs):
        sorted_ccs = sorted(ccs[regType], key = lambda x: twin_age[x['id']])
        if regType.startswith('repeat'):
            repTypeNo += 1
            ax.plot([twin_age[cc['id']] for cc in sorted_ccs],[cc['cc'][0] for cc in sorted_ccs],'o-' if repTypeNo <= 7 else 'o--', label = regType.split('|')[1])

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.set_ylabel('Pearson correlation')
    ax.set_xlabel('Age')

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    fig.savefig(outd('correlation_repeats.png'))

    fig = plt.figure()
    ax = plt.subplot(111)

    regTypes = list(sorted(summary_cc, key = lambda k: summary_cc[k][0]))
#    print regTypes
    _ar = pylab.arange(len(regTypes))
    ax.bar(_ar, [summary_cc[rt][0] for rt in regTypes])
    pylab.xticks( _ar+0.3,  tuple([r.split('|')[-1] if r != 'repeat|All' else r for r in regTypes]), rotation= 'vertical' )
    box = ax.get_position()
    ax.set_position([box.x0, box.y0*3, box.width, box.height*0.7])
    ax.set_ylabel('Pearson correlation')
    ax.set_xlabel('Age')

    fig.savefig(outd('correlation_ranks.png'))


