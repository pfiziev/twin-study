from itertools import izip
import random
from numpy.linalg import *
import numpy as np
from utils import *
import json

if __name__ == '__main__':
    frags = json.load(open(os.path.join(DATA_DIR,'fragments_for_gay_prediction_top100.out')))
    all_frags = set(chain.from_iterable(frags.values()))
    vals = {}
    for fname in os.listdir(DATA_DIR):
        if fname.endswith('.CGmap.regions'):
            tw = fname.split('_')[0]
            vals[tw] = {}
            for l in open(os.path.join(DATA_DIR, fname), 'r'):
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
        for testTwin in vals:
            cfrags = frags[testTwin][:frags_to_take]
            overlap = set(cfrags) if overlap is None else overlap & set(cfrags)

            X = [[vals[tw][fNo] for fNo in cfrags] for tw in sorted(vals) if tw != testTwin]
#            AGE = [1 if tw in GAY_TWINS else 0 for tw in sorted(vals) if tw != testTwin]
            AGE = [twin_random[tw] for tw in sorted(vals) if tw != testTwin]

            x, residues, rank, s = lstsq(X, AGE)

            predAge = np.inner(x, [vals[testTwin][fNo] for fNo in cfrags])
#            print testTwin, predAge, 1 if testTwin in GAY_TWINS else 0
            print testTwin, predAge, twin_random[testTwin]
#            error += abs((1 if predAge >= 0.5 else 0) - (1 if testTwin in GAY_TWINS else 0))
            error += abs( predAge - twin_random[testTwin])
        error = float(error)/ len(twin_age)
        print 'frags',frags_to_take, 'error', error, 'overlap', len(overlap)
        errors.append((error, frags_to_take, len(overlap)))

    elapsed('lstsq')

    print 'min error', min(errors)
