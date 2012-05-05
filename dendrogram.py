from hcluster import *
import sys
sys.path.append('/home/pf/hoffman2/twindata')
from utils import *


def read_data(fname):
    fragments = {}
    for l in open(fname, 'r'):
        buf = l.split('\t')
        fragments[buf[1]] = float(buf[5])
    return fragments

fs = {}

for twin in datafiles:
    fs[twin] = read_data(datafiles[twin] + '.regions')


cut = None
for t in fs:
    cut = cut & set(fs[t].keys()) if cut is not None else set(fs[t].keys())
    print len(cut)
cut = sorted(cut)

M = [[fs[tw][f] for f in cut] for tw in twins]

Y = pdist(M, 'correlation')
Z = average(Y)
dendrogram(Z, labels = tuple(twins), color_threshold=0)
