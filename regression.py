from itertools import izip
import json
import os, glob
import cPickle as pickle
import pprint
import re
import datetime
import sys
from utils import *
from numpy.linalg import *

if __name__ == '__main__':
    DATA_DIR = os.path.join(os.getcwd(), 'data')
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
    X = [[abs(data[t1][s] - data[t2][s]) for s in sites] for t1, t2 in twinpairs]
    AGE = [twin_age[t1] for t1, t2 in twinpairs]

    x, residues, rank, s = lstsq(X,AGE)
    out = open(os.path.join(DATA_DIR, 'regression.out'), 'w')
    out.write('\n'.join('%d\t%f' % (regId, x_i) for x_i, regId in sorted(izip(x, sites))))
    out.close()
    import random

    RAGE = AGE
    random.shuffle(RAGE)
    x, residues, rank, s = lstsq(X,RAGE)
    out = open(os.path.join(DATA_DIR, 'regression_rand.out'), 'w')
    out.write('\n'.join('%d\t%f' % (regId, x_i) for x_i, regId in sorted(izip(x, sites))))
    out.close()


#    json.dump({'X': X, 'AGE': AGE, 'x' : list(x), 'residues': list(residues), 'rank' : rank, 's' : list(s)},
#              open(os.path.join(DATA_DIR, 'regression.out'), 'w'))
    elapsed('lstsq')
