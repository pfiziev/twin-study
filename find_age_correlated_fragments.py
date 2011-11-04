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



getvf = lambda data, s: [data[t][s] for t in TWINS]
getvd = lambda data, s: [abs(data[t1][s] - data[t2][s]) for t1, t2 in twinpairs]



def gen_random_ages():
    random_ages = []
    for i in range(100):
        _x = twin_age.values()
        random.shuffle(_x)
        random_ages.append(_x)
    return random_ages


def cc(sites, data, ages, deltas = False):
    res = {}
    for s in sites:
        values = (getvd if deltas else getvf)(data, s)
        if all(v == 0 for v in values):
            continue

        res[s] = {}
        res[s]['score'], res[s]['pval'] = pearsonr(values, ages)
        res[s]['data'] = json.dumps(values)
    return res

def eucl(vector):
    return math.sqrt(sum(v**2 for v in vector))

def _cosined(values, eucl_values, ages, eucl_ages):
    return sum(v*a for v, a in izip(values, ages))/(eucl_values * eucl_ages)

def cosined(sites, data, ages, deltas = False):
    res = {}
    eucl_ages = eucl(ages)
    random_ages = list(permutations(ages)) if deltas else gen_random_ages()
    for s in sites:
        values = (getvd if deltas else getvf)(data, s)
        if all(v == 0 for v in values):
            continue

        res[s] = {}
        eucl_values = eucl(values)
        res[s]['score'] = _cosined(values, eucl_values, ages, eucl_ages)
        res[s]['pval'] = len(filter(lambda x: x >= res[s]['score'],
                                    [ _cosined(values, eucl_values, rage, eucl_ages) for rage in random_ages])) / float(len(random_ages))
        res[s]['data'] = json.dumps(values)
    return res

def mean(array):
    return float(sum(array))/len(array)

def _dotprod(values, mean_values, ages, mean_ages):
    return sum((val - mean_values) * (age - mean_ages) for val, age in izip(values, ages))


def dotprod(sites, data, ages, deltas = False):
    res = {}
    mean_ages = mean(ages)
    random_ages = list(permutations(ages)) if deltas else gen_random_ages()
    for s in sites:
        values = (getvd if deltas else getvf)(data, s)
        if all(v == 0 for v in values):
            continue
        res[s] = {}
        mean_values = mean(values)
        res[s]['score'] = _dotprod(values, mean_values, ages, mean_ages)
        res[s]['pval'] = len(filter(lambda x: x >= res[s]['score'],
                                    [  _dotprod(values, mean_values, rage, mean_ages) for rage in random_ages])) / float(len(random_ages))
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
    method = cc
    deltas = False
#    method = pearsonr_on_individual_fragments
#    method = pearsonr_on_individual_fragments_minus_means
    res = method(sites, data, TWIN_PAIR_AGES if deltas else TWIN_AGES, deltas = deltas)
    fdr = FDR(res)
#    fdr = 1
    print fdr

    elapsed('calculating scores')

    print 'Accepted:', len([s for s in sorted(res, key = lambda k: res[k]['score'], reverse = True) if res[s].get('pval', -1) <= fdr])

    out = open(os.path.join(DATA_DIR, '%s_%s.out' % (method.func_name, 'deltas' if deltas else 'fragments')), 'w')
    out.write("Fragment number\tscore\tP-value\tData\tAnnotation\n")
    for s in sorted(res, key = lambda k: res[k]['score'], reverse = True):
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



    elapsed('%s_%s.out' % (method.func_name, 'deltas' if deltas else 'fragments'))
