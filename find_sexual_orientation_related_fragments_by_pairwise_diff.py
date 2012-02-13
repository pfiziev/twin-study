from itertools import izip
import json
import os
from pprint import pformat
import random
import re
from utils import *
from scipy import stats


def t_score(group1, group2):
    xbar = mean(group1)
    ybar = mean(group2)
    sx = sum((x - xbar)**2 for x in group1)/(len(group1) - 1)
    sy = sum((y - ybar)**2 for y in group2)/(len(group2) - 1)
    if sx == sy == 0:
        sx = 10**-10
        sy = 10**-10
        if xbar != ybar:
            print'\nAAA'
            print xbar, ybar, sx, sy
            print group1
            print group2


    sp = ((len(group1) - 1)*sx + (len(group2) - 1)*sy)/(len(group1) + len(group2) - 2)
    return (xbar - ybar)/(math.sqrt(sp*(1./len(group1) + 1./len(group2))))


def pw_score(group1, group2):

    diffs = [group1[tw1] - group2[brother[tw1]] for tw1 in group1]

    return sum(diffs) if (abs(max(group1.values()) - min(group2.values()) < 0.01) or abs(max(group2.values()) - min(group1.values())) < 0.01 ) else 0
#    return abs(sum(diffs)) if (all(d > 0 for d in diffs) or all(d < 0 for d in diffs)) else 0

#    return abs(sum(group1[tw1] - group2[brother[tw1]] for tw1 in group1))






def FDR(res, alpha = 0.05):
    pval = -1
    for i, zscore in enumerate(sorted(res, key = lambda zscore: zscore['pval'])):
        if zscore['pval'] <= (i+1)*alpha/len(res):
            pval = zscore['pval']
    return pval


if __name__ == '__main__':
    DATA_DIR = os.path.join(os.getcwd(), 'data')

    # read annotation
    anno = {}
    for l in open(os.path.join(ANNO_DIR,'RRBS_mapable_regions.info.annotated')):
        chrNo, regId, st, en, regAnno = l.split('\t')
        regId = int(regId)
        anno[regId] = {
            'anno'  : json.loads(regAnno),
            'chrNo' : chrNo,
            'fragId' : regId,
            'start' : int(st),
            'end'   : int(en)
        }
    elapsed('annotation')


    # determine common fragments
    common_frags = None
    data = {}
    for fname in os.listdir(DATA_DIR):
        if fname.endswith('.CGmap.fragments_and_sites'):

            cdata = {}

            print fname

            current_frags = set()
            for l in open(os.path.join(DATA_DIR, fname), 'r'):
                chrNo, regId, start, end, regSites = l.split('\t') #filter(None,re.split(r'\t',l))
                regId = int(regId)
                current_frags.add(regId)
                cdata[regId] = json.loads(regSites)

            common_frags = current_frags if common_frags is None else common_frags & current_frags
            data[fname.split('_')[0]] = cdata
            print len(current_frags), len(common_frags)

    common_frags = sorted(common_frags)
    elapsed('reading fragment data')

    mu = 0
    total_sites = 0


    z_score_info = {}
    _c = 0
    _frags = []
    # calculate the t-scores and mu
    for frag in common_frags:
        common_sites = None
        frag_t_scores = []
        meth_levels = dict((tw, dict((s['pos'], s['methLevel']) for s in data[tw][frag])) for tw in data)


        am_gay      = sum(sum(meth_levels[tw][site] for site in meth_levels[tw])/len(meth_levels[tw]) for tw in GAY_TWINS)/len(GAY_TWINS)
        am_hetero   = sum(sum(meth_levels[tw][site] for site in meth_levels[tw])/len(meth_levels[tw]) for tw in HETERO_TWINS)/len(HETERO_TWINS)

        if abs(am_gay - am_hetero) < 0.1:
            continue


        for tw in data:
            common_sites = set(s['pos'] for s in data[tw][frag]) if common_sites is None else common_sites & set(s['pos'] for s in data[tw][frag])

        if len(common_sites) < 2:
            continue

        _c += 1
        if 100 <= _c <= 105:
            _frags.append(frag)
            print '\n' + '+'*100
            print frag
            print pformat(meth_levels)
            print pformat(common_sites)
            print am_gay, am_hetero

        for site in common_sites:
            if 100 <= _c <= 105:
                print site, [meth_levels[tw][site] for tw in GAY_TWINS],[meth_levels[tw][site] for tw in HETERO_TWINS]

            site_t_score = pw_score(dict((tw, meth_levels[tw][site]) for tw in GAY_TWINS),dict((tw, meth_levels[tw][site]) for tw in HETERO_TWINS))
            frag_t_scores.append(site_t_score)
            mu += site_t_score

        # calculate randomized t-bar and sigma
        _mls = {}
        for tw in meth_levels:
            for site in meth_levels[tw]:
                if site not in _mls: _mls[site] = {}
                _mls[site][tw] = meth_levels[tw][site]

        rt_bar_sigma = []
        for rand_i in xrange(100):
            for site in _mls:
                values = _mls[site].values()
                keys = _mls[site].keys()
                random.shuffle(values)
                _mls[site] = dict((tw, val) for tw, val in izip(keys, values))


            for tw in meth_levels:
                for site in meth_levels[tw]:
                    meth_levels[tw][site] = _mls[site][tw]


            rt_scores = [pw_score(dict((tw, meth_levels[tw][site]) for tw in GAY_TWINS),
                                  dict((tw, meth_levels[tw][site]) for tw in HETERO_TWINS)) for site in common_sites]

            if 100 <= _c <= 105 and rand_i == 0:
                print 'RANDOM METH LEVELS'
                print pformat(meth_levels)
                print rt_scores, mean(rt_scores), std(rt_scores)

            rt_bar_sigma.append((mean(rt_scores), std(rt_scores) or 10**-1))


        if not (all(f > 0 for f in frag_t_scores) or all(f < 0 for f in frag_t_scores)):
            continue

        total_sites += len(common_sites)

        z_score_info[frag] = {}
        z_score_info[frag]['tbar'] = mean(frag_t_scores)
        z_score_info[frag]['sigma'] = std(frag_t_scores)
        z_score_info[frag]['n'] = len(common_sites)
        z_score_info[frag]['rt_bar_sigma'] = rt_bar_sigma


    mu /= total_sites
    print 'mu:', mu
    print 'total sites:', total_sites
    elapsed('calculate mu')
    # calculate the z scores
    z_scores = []
    for frag in z_score_info:
        if frag in _frags:
            print '\n\nFINAL'
            print frag
            print pformat(z_score_info[frag])
            print (z_score_info[frag]['tbar'] - mu)*math.sqrt(z_score_info[frag]['n'])/z_score_info[frag]['sigma']
            rz_scores = [(r_tbar - mu)*math.sqrt(z_score_info[frag]['n'])/r_sigma for r_tbar, r_sigma in z_score_info[frag]['rt_bar_sigma']]
            print rz_scores
            print float(len(filter(lambda x: abs(x) >= abs(z_score), rz_scores)))/len(rz_scores)

        z_score = (z_score_info[frag]['tbar'] - mu)*math.sqrt(z_score_info[frag]['n'])/z_score_info[frag]['sigma']
        if abs(z_score) >= 2:
            rz_scores = [(r_tbar - mu)*math.sqrt(z_score_info[frag]['n'])/r_sigma for r_tbar, r_sigma in z_score_info[frag]['rt_bar_sigma']]
            z_scores.append(dict(
                            anno[frag],
                            **{ 'z_score'  :   z_score, # z-score
                                'pval'     :   float(len(filter(lambda x: abs(x) >= abs(z_score), rz_scores)))/len(rz_scores)# p-value

                            }))

    elapsed('z-scores')


    fdr = FDR(z_scores, alpha = 1)
    print 'FDR:', fdr
    accepted = sorted([z for z in z_scores if z['pval'] <= fdr], key = lambda x: abs(x['z_score']), reverse = True)
    print 'Accepted:', len(accepted)
    json.dump(accepted,
              open(os.path.join(DATA_DIR, 'fragments_related_to_sexual_orientation_pw_diff.out'), 'w'),
              indent = 2)

    elapsed('dump')