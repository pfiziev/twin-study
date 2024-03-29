from itertools import izip, permutations
import json
import pprint
import random
import re
import math
from scipy.stats.stats import pearsonr, spearmanr
from utils import *


TWIN_AGES = sorted(twin_age.values())
TWINS =  sorted(twin_age, key = lambda x: twin_age[x])
TWIN_PAIR_AGES = [twin_age[t1] for t1, t2 in twinpairs]



#getvf = lambda data, s: [data[t][s] for t in TWINS]
getvf = lambda data, s: [round (data[t][s],22) for t in TWINS]
getvd = lambda data, s: [round(abs(data[t1][s] - data[t2][s]), 22) for t1, t2 in twinpairs]
#getvd = lambda data, s: [abs(data[t1][s] - data[t2][s]) for t1, t2 in twinpairs]



def gen_random_ages(ages = None):
    random_ages = []
    for i in range(10):
        _x = twin_age.values() if ages is None else ages
        random.shuffle(_x)
        random_ages.append(_x)
    return random_ages

def _pearsonr(values, mean_values, ages, mean_ages):
    val_diffs = [v - mean_values for v in values]
    age_diffs = [age - mean_ages for age in ages]
    return float(sum(vald*aged for vald, aged in izip(val_diffs, age_diffs)))/math.sqrt(sum(vald*vald for vald in val_diffs)*sum(aged*aged for aged in age_diffs))



def _mono(values):
    diff = []
    for i in xrange(len(values)-1):
        diff.append(values[i+1] - values[i])
    if all(d > 0 for d in diff) or all(d < 0 for d in diff):
        return sum(diff)
    else:
        return 0


def mono(sites, data, ages, deltas = False):
    res = {}
    random_ages = list(permutations(ages)) if deltas else gen_random_ages()
    mean_ages = mean(ages)
    for s in sites:
        values = (getvd if deltas else getvf)(data, s)
        if len(set(values)) == 1:
            continue

        res[s] = {}
        res[s]['score'] =  _mono(values)

        res[s]['pval'] = 1
#        res[s]['score'], res[s]['pval'] = pearsonr(values, ages)




#        mean_values = mean(values)
#        res[s]['score'], res[s]['pval'] = _pearsonr(values, mean_values, ages, mean_ages), 0
        res[s]['data'] = json.dumps(values)

#        res[s]['rscore'] = json.dumps([ _pearsonr(values, mean_values, rage, mean_ages) for rage in random_ages])
#        rscores = [ pearsonr(values,rage) for rage in random_ages]
        res[s]['rscore'] = [_mono(sorted(values))] #json.dumps([r[0] for r in rscores])
        res[s]['rpvals'] = '' #json.dumps([r[1] for r in rscores])
    return res

def cc(sites, data, ages, deltas = False):
    res = {}
#    random_ages = list(permutations(ages)) if deltas else gen_random_ages()
    random_ages =  gen_random_ages(ages)
#    mean_ages = mean(ages)
    for s in sites:
        values = (getvd if deltas else getvf)(data, s)
        if len(set(values)) == 1:
            continue

        if max(values) - min(values) < 0.1:
            continue
        res[s] = {}

        res[s]['score'], res[s]['pval'] = pearsonr(values, ages)
#        mean_values = mean(values)
#        res[s]['score'], res[s]['pval'] = _pearsonr(values, mean_values, ages, mean_ages), 0
        res[s]['data'] = json.dumps(values)

#        res[s]['rscore'] = json.dumps([ _pearsonr(values, mean_values, rage, mean_ages) for rage in random_ages])
        rscores = [ pearsonr(values,rage) for rage in random_ages]
        res[s]['rscore'] = json.dumps([r[0] for r in rscores])
        res[s]['rpvals'] =  json.dumps([r[1] for r in rscores])
    return res


def spearmanc(sites, data, ages, deltas = False):
    res = {}
    random_ages = list(permutations(ages)) if deltas else gen_random_ages()

    for s in sites:
        values = (getvd if deltas else getvf)(data, s)
        if len(set(values)) == 1:
            continue

        res[s] = {}
        res[s]['score'], res[s]['pval'] = spearmanr(values, ages)

        res[s]['rscore'] = json.dumps([spearmanr(values, rage)[0] for rage in random_ages])


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

        rscores = [ _cosined(values, eucl_values, rage, eucl_ages) for rage in random_ages]

        res[s]['pval'] = len(filter(lambda x: x >= res[s]['score'], rscores)) / float(len(random_ages))

#        res[s]['rscore'] = sum(rscores)/ float(len(random_ages))
        res[s]['rscore'] = json.dumps(rscores)

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
    rpvals = [float(i)/120 for i in xrange(1,121)]
    for s in sites:
        values = (getvd if deltas else getvf)(data, s)
        if all(v == 0 for v in values):
            continue
        res[s] = {}
        mean_values = mean(values)
        res[s]['score'] = _dotprod(values, mean_values, ages, mean_ages)

        rscores = [  _dotprod(values, mean_values, rage, mean_ages) for rage in random_ages]

        res[s]['pval'] = len(filter(lambda x: x >= res[s]['score'], rscores)) / float(len(random_ages))

#        res[s]['rscore'] = sum(rscores)/ float(len(random_ages))

        res[s]['rscore'] = '' #json.dumps(rscores)

        res[s]['rpvals'] = '' # json.dumps(rpvals)
        res[s]['data'] = json.dumps(values)
    return res


def FDR(res, alpha = 0.05):
    pval = -1
    for i, s in enumerate(sorted(res, key = lambda k: res[k]['pval'])):
#        if i < 10:
#            print i,s, pprint.pformat(res[s])
        if res[s]['pval'] <= (i+1)*alpha/len(res):
            pval = res[s]['pval']
#        else:
#            break
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
    for twin_id in datafiles:
        fname = reg_fname(twin_id)
        cdata = {}

        print fname

        current_sites = set()
        for l in open(fname):
            chrNo, regId, start, end, regSites, methLevel = filter(None,re.split(r'\s+',l))
            regId = int(regId)
            current_sites.add(regId)
            cdata[regId] = float(methLevel)

        sites = current_sites if sites is None else sites & current_sites
        data[twin_id] = cdata
        print len(current_sites), len(sites)

    elapsed('reading data')
    sites = sorted(sites)

#    res = pearsonr_on_individual_fragment_deltas(sites, data)
#    method = dotprod_on_individual_fragments
    method = cc
    deltas = True
    alpha = 0.15
#    method = pearsonr_on_individual_fragments
#    method = pearsonr_on_individual_fragments_minus_means
    res = method(sites, data, TWIN_PAIR_AGES if deltas else TWIN_AGES, deltas = deltas)
#    fdr = FDR(res, alpha=0.37)
    fdr = FDR(res, alpha = alpha)
    fdr = 1
    print fdr

    elapsed('calculating scores')

    print 'Accepted:', len([s for s in sorted(res, key = lambda k: res[k]['score'], reverse = True) if res[s].get('pval', -1) <= fdr])

    out = open(outd('%s_%s_%.2lf.out' % ('deltas' if deltas else 'fragments', method.func_name, alpha)), 'w')
    out.write("#fragment\tscore\tp-value\trscore\tdata\tannotation\n")
    for s in sorted(res, key = lambda k: res[k]['score'], reverse = True):
        if res[s].get('pval', -1) <= fdr:
            out.write("%d\t%f\t%f\t%s\t%s\t%s\t%s\n" % (s, res[s]['score'], res[s].get('pval', -1), res[s]['rscore'],res[s].get('rpvals',''), res[s]['data'], anno[s]))
    out.close()

    elapsed('real data')

#    rtimes = 1
#    raccepted = 0
#    rage = list(reversed(TWIN_PAIR_AGES if deltas else TWIN_AGES))
#
#    for rs in xrange(rtimes):
##        random.shuffle(rage)
#        res = method(sites, data, rage, deltas = deltas)
#        fdr = 1
#
#
#        out = open(os.path.join(DATA_DIR, '%s_%s_%.2lf.random_out' % ('deltas' if deltas else 'fragments', method.func_name, alpha)), 'w')
#        out.write("#fragment\tscore\tp-value\tdata\tannotation\n")
#
#        for s in sorted(res, key = lambda k: res[k]['score'], reverse = True):
#            if res[s].get('pval', -1) <= fdr:
#                out.write("%d\t%f\t%f\t%s\t%s\n" % (s, res[s]['score'], res[s].get('pval', -1), res[s]['data'], anno[s]))
#        out.close()
#


#        fdr = FDR(res, alpha = alpha)
#        print pprint.pformat([s for s in sorted(res, key = lambda k: res[k]['score'], reverse = True) if res[s].get('pval', -1) <= fdr][:10])
#        print 'random fdr:', fdr
#        _raccepted =  len([s for s in sorted(res, key = lambda k: res[k]['score'], reverse = True) if res[s].get('pval', -1) <= fdr])
#        print '_raccepted:', _raccepted, '\n\n'
#        raccepted += _raccepted
#    print "Random: ", float(raccepted)/rtimes






    elapsed('%s_%s.out' % (method.func_name, 'deltas' if deltas else 'fragments'))
