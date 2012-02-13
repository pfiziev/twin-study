import json
import os
from utils import *

__author__ = 'pf'


if __name__ == '__main__':
    top = 3751
#    top = 375000
    lno = 0
    genes = {}
    for l in open(os.path.join(DATA_DIR, 'deltas_cc_1.00.out')):
        if l[0] == '#': continue

        fno, score, pval, rscore, rpval, data, anno = l.split('\t')

        for a in json.loads(anno):
            if a['type'] in ['gene', 'promotor']:
                if a['info'] not in genes: genes[a['info']] = []
                genes[a['info']].append({'fno' : fno, 'score' : score, 'data' : data})

        lno += 1
        if lno >= top:
            print 'Min score: ', score, ' pval:',pval
            break

    cnt = 0
    out = open(os.path.join(DATA_DIR,'age_related_genes.out'), 'w')
    for g in filter(lambda k: len(genes[k]) >= 2, reversed(sorted(genes, key = lambda k: len(genes[k])))):
        cnt += 1
#        print '%s\t%d\t%s' % (g, len(genes[g]), ','.join(a['fno'] for a in genes[g]))
        print '%s\t%d\t%s' % (g, len(genes[g]), ','.join('[%s:%s]' % (a['fno'], a['score']) for a in genes[g]))
        out.write('%s\n'%g)

    out.close()
    print cnt