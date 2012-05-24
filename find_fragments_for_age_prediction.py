from itertools import izip, permutations
import json
import multiprocessing
import pprint
import random
import re
import math
from scipy.stats.stats import pearsonr
from utils import *



def _pearsonr(val_diffs, age_diffs):
    return float(sum(vald*aged for vald, aged in izip(val_diffs, age_diffs))) #/math.sqrt(sum(vald*vald for vald in val_diffs)*sum(aged*aged for aged in age_diffs))


TWINS =  sorted(twin_age, key = lambda x: twin_age[x])

def _cc(sites, data, ages, testTwin = None):
    res = {}
    mean_ages = mean(ages)

    random_agediffs = []

    for i in range(10000):
        _x = list(ages)
        random.shuffle(_x)
        random_agediffs.append([age - mean_ages for age in _x])

    age_diffs = [age - mean_ages for age in ages]

    for s in sites:
        values = [round(data[t][s],22) for t in TWINS if t != testTwin]
        if len(set(values)) == 1:
            continue

        if max(values) - min(values) < 0.1:
            continue
        res[s] = {}

#        res[s]['score'], res[s]['pval'] = pearsonr(values, ages)
        mean_values = mean(values)
        val_diffs = [v - mean_values for v in values]
        res[s]['score'] = _pearsonr(val_diffs, age_diffs)
        res[s]['pval'] = len(filter(lambda x: abs(x) >= abs(res[s]['score']),
                                    [ _pearsonr(val_diffs, rage_diffs) for rage_diffs in random_agediffs])) / float(len(random_agediffs))



        res[s]['data'] = json.dumps(values)

#        rscores = [ pearsonr(values,rage) for rage in random_ages]
#        res[s]['rscore'] = json.dumps([r[0] for r in rscores])
#        res[s]['rpvals'] =  json.dumps([r[1] for r in rscores])
    return res

def cc(sites, data, ages, testTwin = None):
    res = {}

    for s in sites:
        values = [round(data[t][s],22) for t in TWINS if t not in [testTwin, brother[testTwin]]]
        
        if len(set(values)) == 1:
            continue

        if max(values) - min(values) < 0.1:
            continue
        res[s] = {}

        res[s]['score'], res[s]['pval'] = pearsonr(values, ages)
   
        res[s]['data'] = json.dumps(values)

    return res




def FDR(res, alpha = 0.05):
    pval = -1
    for i, s in enumerate(sorted(res, key = lambda k: res[k]['pval'])):
        if res[s]['pval'] <= (i+1)*alpha/len(res):
            pval = res[s]['pval']
    return pval

#
#def process((sites, data, testTwin)):
#    print 'starting', testTwin
#    highpriority()
#    DATA_DIR = os.path.join(os.getcwd(), 'data')
#    alpha = 0.05
#    res = cc(sites, data, [twin_age[tw] for tw in TWINS if tw != testTwin], testTwin = testTwin)
#
#    fdr = FDR(res, alpha = alpha)
#
#    result = [s for s in sorted(res, key = lambda k: res[k]['score'], reverse = True) if res[s]['pval'] <= fdr]
#
#    print testTwin, 'FDR:', fdr, 'Accepted:', len(result)
#
#    json.dump(result, open(os.path.join(DATA_DIR, '%s_%.2lf_fragments_for_age_prediction.out' % (testTwin, alpha)), 'w'))
#
#    elapsed(testTwin)


if __name__ == '__main__':

    sites = None
    data = {}
    for twin_id, fname in datafiles.items():
        fname += '.regions'
        cdata = {}
        print fname

        current_sites = set()

        for l in open(fname, 'r'):
            chrNo, regId, start, end, regSites, methLevel = filter(None, re.split(r'\s+',l))
            regId = int(regId)
            current_sites.add(regId)
            cdata[regId] = float(methLevel)

        sites = current_sites if sites is None else sites & current_sites
        data[twin_id] = cdata
        print len(current_sites), len(sites)

    elapsed('reading data')
    sites = sorted(sites)

    alpha = 0.05

#    multiprocessing.Pool(processes = 6).map(process, [(sites, data, testTwin) for testTwin in data])


    result = {}
    for testTwin in sorted(data):

        res = cc(sites, data, [twin_age[tw] for tw in TWINS if tw not in [testTwin, brother[testTwin]]], testTwin = testTwin)

#        fdr = FDR(res, alpha = alpha)

#        result[testTwin] = [s for s in sorted(res, key = lambda k: res[k]['score'], reverse = True) if res[s].get('pval', -1) <= fdr]

        result[testTwin] = [(s, res[s]['score']) for s in sorted(res, key = lambda k: abs(res[k]['score']), reverse = True)][:100]

#        print testTwin , 'FDR:', fdr, 'Accepted:', len(result[testTwin])

        elapsed(testTwin)

    json.dump(result, open(os.path.join(DATA_DIR, 'fragments_for_age_prediction_top100.out_new'), 'w'))

    elapsed('end')
