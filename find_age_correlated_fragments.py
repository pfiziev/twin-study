from itertools import izip, permutations
import json
import pprint
import random
import re
import math
from scipy.stats.stats import pearsonr
from utils import *


TWIN_AGES = sorted(twin_age.values())
TWINS =  sorted(twin_age, key = lambda x: twin_age[x])
TWIN_PAIR_AGES = [twin_age[t1] for t1, t2 in twinpairs]


def mean(array):
    return float(sum(array))/len(array)

def dotprod(v1, ages, eucl_ages = None, mean_ages = None):
    m1 = mean(v1)
    return sum((val1 - m1) * (val2 - mean_ages)  for val1, val2 in izip(v1, ages))

def eucl(vector):
    return math.sqrt(sum(v**2 for v in vector))

def cosined(v1, ages, eucl_ages = None, mean_ages = None, eucl_v1 = None):
    return sum(v*a for v, a in izip(v1, ages))/((eucl_v1 if eucl_v1 else eucl(v1)) * eucl_ages)

def calc_pval(method, score, values, random_ages = None):
    if not hasattr(calc_pval, 'RTWIN_AGES'):
        if random_ages is None:
            random_ages = []
            for i in range(100):
                _x = twin_age.values()
                random.shuffle(_x)
                random_ages.append(_x)

        calc_pval.RTWIN_AGES = random_ages
    
    return len(filter(lambda x: x >= score,
                [ method(values, rage) for rage in calc_pval.RTWIN_AGES])) / float(len(calc_pval.RTWIN_AGES))



def apply_method(method, values, ages, random_ages = None):
    if not hasattr(apply_method, 'cache'):
        apply_method.cache = {'eucl_ages' : eucl(ages),
                              'mean_ages' : mean(ages)}
    if not all(v == 0 for v in values):
        if method is pearsonr:
            score, pval = pearsonr(values, ages)
        else:
            score = method(values, ages, **apply_method.cache)
            pval = calc_pval(method, score, values, random_ages = random_ages)
    else:
        score = 0
        pval = 1
    return score, pval


def deltas(method, sites, data, ages = None):
    if ages is None: ages = TWIN_PAIR_AGES

    random_ages = list(permutations(ages))
    res = {}
    for s in sites:
        deltas = [abs(data[t1][s] - data[t2][s]) for t1, t2 in twinpairs]
        res[s] = {}
        res[s]['score'], res[s]['pval'] = apply_method(method, deltas, ages, random_ages = random_ages)

        res[s]['data'] = json.dumps({'deltas': deltas, 'data' : [[data[t1][s], data[t2][s]] for t1, t2 in twinpairs]})
    return res



def fragments(method, sites, data, ages = None):
    if ages is None: ages = TWIN_AGES

    res = {}

    for s in sites:
        values = [data[t][s] for t in TWINS]
        res[s] = {}
        res[s]['score'], res[s]['pval'] = apply_method(method, values, ages)
        res[s]['data'] = json.dumps(values)
    return res



def FDR(res, alpha = 0.05):
    pval = -1
    for i, s in enumerate(sorted(res, key = lambda k: res[k]['pval'])):
#        if i < 10:
#            print i,s, pprint.pformat(res[s])
        if res[s]['pval'] <= (i+1)*alpha/len(res):
            pval = res[s]['pval']
        else:
            break
    return pval
        

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
#    res = pearsonr_on_individual_fragment_deltas(sites, data)
#    method = dotprod_on_individual_fragments
    resolution = fragments
    method = cosined
#    method = pearsonr_on_individual_fragments
#    method = pearsonr_on_individual_fragments_minus_means
    res = resolution(method, sites, data)
    fdr = FDR(res)
#    fdr = 1
    print fdr

    elapsed('calculating scores')

    print 'Accepted:', len([s for s in sorted(res, key = lambda k: res[k]['score'], reverse = True) if res[s].get('pval', -1) <= fdr])

    out = open(os.path.join(DATA_DIR, '%s_%s.out' % (resolution.func_name, method.func_name)), 'w')
    out.write("Fragment number\tscore\tP-value\tData\tAnnotation\n")
    for s in sorted(res, key = lambda k: res[k]['score'], reverse = True):
        if s == 147302:
            print "%d\t%f\t%f\t%s\t%s\n" % (s, res[s]['score'], res[s].get('pval', -1), res[s]['data'], anno[s])
        if res[s].get('pval', -1) <= fdr:
            out.write("%d\t%f\t%f\t%s\t%s\n" % (s, res[s]['score'], res[s].get('pval', -1), res[s]['data'], anno[s]))
    out.close()

    elapsed('real data')
    
#    rtimes = 2
#    raccepted = 0
#    rage = TWIN_AGES
#
#    for rs in xrange(rtimes):
#        random.shuffle(rage)
#        res = resolution(method, sites, data, ages = rage)
#        fdr = FDR(res)
#        print pprint.pformat([s for s in sorted(res, key = lambda k: res[k]['score'], reverse = True) if res[s].get('pval', -1) <= fdr][:10])
#        print 'random fdr:', fdr
#        _raccepted =  len([s for s in sorted(res, key = lambda k: res[k]['score'], reverse = True) if res[s].get('pval', -1) <= fdr])
#        print '_raccepted:', _raccepted, '\n\n'
#        raccepted += _raccepted
#    print "Random: ", float(raccepted)/rtimes
#



    elapsed('%s_%s.out' % (resolution.func_name, method.func_name))
