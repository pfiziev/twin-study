from itertools import izip
import random
from numpy.linalg import *
import numpy as np
from utils import *
import json

if __name__ == '__main__':
    frags = json.load(open(os.path.join(DATA_DIR,'fragments_for_age_prediction_top100.out')))
    all_frags = set(f[0] for fvals in frags.values() for f in fvals)
    vals = {}
    for tw in datafiles:
        fname = reg_fname(tw)
        vals[tw] = {}
        for l in open(fname, 'r'):
            _, fNo, _, _, _, methLevel = l.split('\t')
            fNo = int(fNo)
            if fNo in all_frags:
                vals[tw][fNo] = float(methLevel)
        elapsed(fname)

    errors = []
    overlaps = []
    for frags_to_take in range(1, 1 + len(frags.values()[0])):
#    for frags_to_take in range(37, 38):
        error = 0
        overlap = None
        for testTwin in sorted(vals):
            cfrags = [f[0] for f in frags[testTwin][:frags_to_take]]
            overlap = set(cfrags) if overlap is None else overlap & set(cfrags)

            X = [[vals[tw][fNo] for fNo in cfrags] for tw in sorted(vals) if tw not in [testTwin, brother[testTwin]]]
            AGE = [twin_age[tw] for tw in sorted(vals) if tw not in [testTwin, brother[testTwin]]]
#
#            X = [[vals[tw][fNo] for fNo in cfrags] for tw in sorted(vals) if tw not in [testTwin]]
#            AGE = [twin_age[tw] for tw in sorted(vals) if tw not in [testTwin]]

            x, residues, rank, s = lstsq(X, AGE)

            predAge = np.inner(x, [vals[testTwin][fNo] for fNo in cfrags])
            print testTwin, predAge, twin_age[testTwin]
            error += (predAge - twin_age[testTwin])**2
        error = math.sqrt(float(error) / len(twin_age))
        print 'frags',frags_to_take, 'error', error, 'overlap', len(overlap)
        errors.append((error, frags_to_take, set(overlap)))

    rmse, frags_to_take, overlap = min(errors)
    # Build a linear model from the fragments that are common to all twins in the minumum
    error = 0
    for testTwin in sorted(vals):

        X = [[vals[tw][fNo] for fNo in overlap] for tw in sorted(vals) if tw not in [testTwin, brother[testTwin]]]
        AGE = [twin_age[tw] for tw in sorted(vals) if tw not in [testTwin, brother[testTwin]]]

        x, residues, rank, s = lstsq(X, AGE)

        predAge = np.inner(x, [vals[testTwin][fNo] for fNo in overlap])
        print testTwin, predAge, twin_age[testTwin]
        error += (predAge - twin_age[testTwin])**2
    error = math.sqrt(float(error) / len(twin_age))
    print 'final min error:', error




    elapsed('lstsq')
    open(outd('regression.out'), 'w').write('min RMSE:' + str(rmse)+' years\tfragments:'+str(frags_to_take)+'\toverlap:'+str(len(overlap))+'\n')
    print 'min error', min(errors)
