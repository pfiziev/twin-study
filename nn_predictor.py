import os
from pybrain.tools.shortcuts import buildNetwork
from pybrain.datasets import SupervisedDataSet
from pybrain.supervised.trainers import BackpropTrainer
import random
from utils import *

if __name__ == '__main__':
    frags = set()
    for l in open(os.path.join(DATA_DIR,'fragments_cc_0.05.out')):
        if l.startswith('#'): continue
        frags.add(l.split('\t')[0])
    vals = {}
    for fname in os.listdir(DATA_DIR):
        if fname.endswith('.CGmap.regions'):
            tw = fname.split('_')[0]
            vals[tw] = {}
            for l in open(os.path.join(DATA_DIR, fname), 'r'):
                _, fNo, _, _, _, methLevel = l.split('\t')
                if fNo in frags:
                    vals[tw][fNo] = float(methLevel)
            elapsed(fname)

    input_size = len(frags)
    print 'input size', input_size
    net = buildNetwork(input_size, input_size, 1)
    elapsed('building nn')
    ds = SupervisedDataSet(input_size, 1)
    testTwin = random.choice(vals.keys())
    for tw in vals:
        if tw != testTwin:
            ds.addSample(tuple(vals[tw][fNo] for fNo in sorted(vals[tw].iterkeys())), (twin_age[tw],))
    elapsed('preparing sample')
    trainer = BackpropTrainer(net, ds)
    tres = trainer.trainUntilConvergence()
    elapsed('training')
    print 'Test twin', testTwin
    print net.activate([vals[testTwin][fNo] for fNo in sorted(vals[testTwin].iterkeys())])
    elapsed('testing')