import json
import re
from scipy.stats.stats import pearsonr
from utils import *

pred_sites = {'cg09809672' : False, 'cg27210390': False, 'cg12799895' : False}

if __name__ == '__main__':
    cgs = {}
    for l in open(os.path.join(DATA_DIR,'plos_paper\journal.pone.0014821.s002.json')):
        row = json.loads(l)
        row['mvals'] = {}
        cgs[tuple(row['position'])] = row
    out = open(os.path.join(DATA_DIR,'plos_paper\journal.pone.0014821.s002.eval_sites'), 'w')
    for fname in os.listdir(DATA_DIR):
        if fname.endswith('.CGmap'):
            cnt = 0
            tw = fname.split('_')[0]
            for l in open(os.path.join(DATA_DIR, fname), 'r'):
                chrNo, nucl, pos, methType, methSubtype, methLevel, methReads, totalReads = filter(None,re.split(r'\s+',l))
                chrNo = int(chrNo)
                chrNo = chrNo if chrNo <= 22 else ['X','Y','M'][chrNo-23]
                key = (chrNo, int(pos))
                if key in cgs:
                    out.write('#' + l)
                    if int(totalReads) >= 4 and methType == 'CG':
                        cgs[key]['mvals'][tw] = float(methLevel)
                        cnt += 1
            elapsed(fname+': ' + str(cnt))
            out.write('#'+fname+': ' + str(cnt)+'\n')
    for k in cgs:
        tws = [tw for tw in twins if tw in cgs[k]['mvals']]
        mvals = [cgs[k]['mvals'][tw] for tw in tws]
        cgs[k]['tw_cc'] = pearsonr(mvals,
                                [twin_age[tw] for tw in tws]) if len(tws) > 1 and len(set(mvals)) > 1 else (0, 1)
    out.write(json.dumps([cgs[k] for k in cgs])+ '\n')
    out.close()

    """
# the code to plot the comparison from dreampie
# plos_sites is loaded from the output of this script

import numpy as np
import matplotlib.pyplot as plt
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
plt.xlabel('CG site (24 sites in total)')

plt.title('Pearson correlation between methylation levels and age')

plt.legend( (rects1[0], rects2[0]), ('This data', 'Vilain et al., PLoS 2011'), loc = 'middle' )

    """