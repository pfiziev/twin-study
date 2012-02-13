import json
import re
from scipy.stats.stats import pearsonr
from utils import *



if __name__ == '__main__':

    cgs = {}
    for l in open(os.path.join(DATA_DIR,'plos_paper\journal.pone.0014821.s002.json')):
        row = json.loads(l)
        row['mvals'] = {}
        if row['fragment'] > 0:
            cgs[row['fragment']] = row


    out = open(os.path.join(DATA_DIR,'plos_paper\journal.pone.0014821.s002.eval_fragments'), 'w')
    for fname in os.listdir(DATA_DIR):
        if fname.endswith('.CGmap.regions'):
            cnt = 0
            tw = fname.split('_')[0]
            for l in open(os.path.join(DATA_DIR, fname), 'r'):
                chrNo, fragId, st, en, reads, methLevel = filter(None,re.split(r'\s+',l))
                fragId = int(fragId)
                if fragId in cgs:
                    out.write('#' + l)
                    cgs[fragId]['mvals'][tw] = float(methLevel)
                    cnt += 1
            elapsed(fname+': ' + str(cnt))
            out.write('#'+fname+': ' + str(cnt)+'\n')

    to_del = []
    for k in cgs:

        tws = [tw for tw in twins if tw in cgs[k]['mvals']]
        if len(tws) == 0:
            to_del.append(k)

        mvals = [cgs[k]['mvals'][tw] for tw in tws]
        cgs[k]['tw_cc'] = pearsonr(mvals, [twin_age[tw] for tw in tws]) if len(tws) > 1 and len(set(mvals)) > 1 else (0, 1)

        print k, cgs[k]['Symbol'], cgs[k]['tw_cc'], cgs[k]['r'] ,len(tws)

    for k in to_del:
        del cgs[k]
    print 'total', len(cgs)
    plos_sites = [cgs[k] for k in cgs]
    out.write(json.dumps([cgs[k] for k in cgs])+ '\n')
    out.close()


    # the code to plot the comparison from dreampie

    import numpy as np
    import matplotlib.pyplot as plt
    fig = plt.figure()

    site_ccs = sorted([[p['tw_cc'][0], p['r']] for p in plos_sites if p['tw_cc'] != [0,1]], key = lambda x: abs(x[0] - x[1]))

    ind = np.arange(len(site_ccs))  # the x locations for the groups
    width = 0.35       # the width of the bars


    plt.subplot(111)
    rects1 = plt.bar(ind, [cc[0] for cc in site_ccs], width,
                        color='r')

    rects2 = plt.bar(ind+width, [cc[1] for cc in site_ccs], width,
                        color='y')

    # add some
    plt.ylabel('Pearson correlation')
    plt.xlabel('Fragment (%d in total)' % len(cgs))

    plt.title('Pearson correlation between methylation levels and age')

    plt.legend( (rects1[0], rects2[0]), ('This data', 'Vilain et al., PLoS 2011'), loc = 'best' )

    fig.savefig('plos_fragments.png')




#    cgs = {}
#    for l in open(os.path.join(DATA_DIR,'plos_paper\journal.pone.0014821.s002.json')):
#        row = json.loads(l)
#        cgs[row['Symbol']] = '+' if row['r'] > 0 else '-'
#
#    gene_cc = {}
#    for l in open(os.path.join(DATA_DIR,'fragments_cc_1.00.out')):
#        if l.startswith('#'): continue
#        buf = l.split('\t')
#        cc = float(buf[1])
#        for anno in json.loads(buf[-1]):
#            if 'geneSymbol' in anno.get('info',{}):
#                key = anno['info']['geneSymbol']
#                if key not in gene_cc:
#                    gene_cc[key] = []
#                if abs(cc) >= 0.5:
#                    gene_cc[key].append(['+' if cc > 0 else '-', cc])
#
#    nf, f, cons = 0, 0, 0
#    for g in cgs:
#
#        if g not in gene_cc: nf += 1
#        else:
#            print g, cgs[g], gene_cc[g],  cgs[g] in [cc[0] for cc in gene_cc[g]]
#            f += 1
#            if cgs[g] in [cc[0] for cc in gene_cc[g]]:
#                cons += 1
#
#    print 'Not found:', nf
#    print 'Found:', f
#    print 'Consistent', cons