import json
import os
import re
from utils import *
from scipy import stats




if __name__ == '__main__':
    DATA_DIR = os.path.join(os.getcwd(), 'data')
    anno = {}
    for l in open(os.path.join(ANNO_DIR,'RRBS_mapable_regions.info.annotated')):
        _, regId, _, _, regAnno = l.split('\t')
        anno[int(regId)] = regAnno.strip()
    elapsed('annotation')

    sites = None
    data = {}
    for fname in os.listdir(DATA_DIR):
        if fname.endswith('.CGmap.regions'):

            cdata = {}

            print fname

            current_sites = set()
            for l in open(os.path.join(DATA_DIR, fname), 'r'):
                chrNo, regId, start, end, regSites, methLevel = filter(None,re.split(r'\s+',l))
                regId = int(regId)
                current_sites.add(regId)
                cdata[regId] = float(methLevel)

            sites = current_sites if sites is None else sites & current_sites
            data[fname.split('_')[0]] = cdata
            print len(current_sites), len(sites)

    elapsed('reading data')
    sites = sorted(sites)

    res = {}
    are_dense = lambda array: max(array) - min(array) <= 0.1
    last_site = None
    for s in sites:
        gay = [round (data[t][s],22) for t in GAY_TWINS]
        hetero = [round (data[t][s],22) for t in HETERO_TWINS]
#        if not (are_dense(gay) and are_dense(hetero)):
#           continue

        if abs(mean(gay) - mean(hetero)) < 0.1:
            continue

        last_site = s
        res[s] = {}
        res[s]['z'] = '%f %f %f' %(mean(gay), mean(hetero), abs(mean(gay) - mean(hetero)))
        res[s]['score'], res[s]['pval'] = abs(mean(gay) - mean(hetero)), 0
#        res[s]['score'],  res[s]['pval'] = stats.ttest_ind(gay, hetero)
        #        res[s]['score'], res[s]['pval'] = _pearsonr(values, mean_values, ages, mean_ages), 0
        res[s]['data'] = json.dumps({'gay' : gay, 'hetero' : hetero})
        res[s]['anno'] = anno[s]
#        res[s]['rscore'] = json.dumps([ _pearsonr(values, mean_values, rage, mean_ages) for rage in random_ages])
#       rscores = [ pearsonr(values,rage) for rage in random_ages]
#       res[s]['rscore'] = json.dumps([r[0] for r in rscores])
#       res[s]['rpvals'] =  json.dumps([r[1] for r in rscores])

    elapsed('T tests')
    out = open(os.path.join(DATA_DIR, 'sexual_orientation_fragments.out'), 'w')
    out.write("#site\t%s\n" % '\t'.join(reversed(sorted(res[last_site].keys()))))

    for s in sorted(res, key = lambda k: res[k]['score']):
       if res[s].get('pval', -1) <= 1:
           out.write("%d\t%s\n" % (s, '\t'.join(str(res[s][k]) for k in reversed(sorted(res[last_site].keys())))))
    out.close()

    elapsed('real data')
