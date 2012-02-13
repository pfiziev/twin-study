from itertools import product, chain
import os
import datetime
import math


_twins = dict(  T1 = 'XB002L1',
                T2 = 'XB002L2',
                T3 = 'XB002L3',
                T4 = 'XB002L5',
                T5 = 'XB002L6',
                T6 = 'VB015L1',
                T7 = 'VB015L2',
                T8 = 'VB015L3',
                T9 = 'XA004L5',
                T10 = 'XA004L6' )


twinpairs = [('T7','T8'),
             ('T5','T6'),
             ('T1','T2'),
             ('T9','T10'),
             ('T3','T4')]

# True if gay, False if hetero
GAY_TWINS    = frozenset(['T1', 'T3', 'T5', 'T7', 'T9' ])
HETERO_TWINS = frozenset(['T2', 'T4', 'T6', 'T8', 'T10'])


random_twin_pairs = [pair for pair in product([tp[0] for tp in twinpairs],[tp[1] for tp in twinpairs]) if pair not in twinpairs]

#brother = {'T1': 'T5', # random brothers
# 'T10': 'T9',
# 'T2': 'T4',
# 'T3': 'T1',
# 'T4': 'T2',
# 'T5': 'T10',
# 'T6': 'T8',
# 'T7': 'T3',
# 'T8': 'T7',
# 'T9': 'T6'}


brother = dict((t1, t2) for t1, t2 in twinpairs + [(t2, t1) for t1, t2 in twinpairs])

twin_age = {'T7'  : 21,
            'T8'  : 21,
            'T5'  : 30,
            'T6'  : 30,
            'T1'  : 38,
            'T2'  : 38,
            'T9'  : 45,
            'T10' : 45,
            'T3'  : 55,
            'T4'  : 55}


twin_random = {'T1': 23,
 'T10': 33,
 'T2': 37,
 'T3': 31,
 'T4': 55,
 'T5': 38,
 'T6': 36,
 'T7': 29,
 'T8': 49,
 'T9': 47}



twins = sorted(twin_age, key = lambda k: twin_age[k])



DATA_DIR = os.path.join(os.getcwd(), 'data') 

_datafiles = [os.path.join(DATA_DIR, fname) for fname in os.listdir(DATA_DIR) if fname.endswith('.CGmap')]
datafiles = dict((key, filter(lambda x: _twins[key] in x, _datafiles)[0]) for key in _twins)

stime = datetime.datetime.now()
def elapsed(msg = None):
    print "[%s]" % msg if msg is not None else "+", "Elapsed: " , datetime.datetime.now() - stime


ANNO_DIR = os.path.join(DATA_DIR,'annotation')

regionsinfo = os.path.join(ANNO_DIR, "RRBS_mapable_regions.info")


def highpriority():
    """ Set the priority of the process to high in windows."""

    import sys
    try:
        sys.getwindowsversion()
    except:
        isWindows = False
    else:
        isWindows = True

    if isWindows:
        # Based on:
        #   "Recipe 496767: Set Process Priority In Windows" on ActiveState
        #   http://code.activestate.com/recipes/496767/
        import win32api, win32process, win32con

        pid = win32api.GetCurrentProcessId()
        handle = win32api.OpenProcess(win32con.PROCESS_ALL_ACCESS, True, pid)
        win32process.SetPriorityClass(handle, win32process.HIGH_PRIORITY_CLASS)


def addTo(_dict, _key, _val, default = None):
    if _key not in _dict:
        _dict[_key] = default if default is not None else 0
    _dict[_key] += _val


def overlap(top = 1000, *args):
    alls = None
    for fname in args:
        lno = 0
        ids = []
        for l in open(fname,'r'):
            if lno >= top: break
            ids.append(l.split('\t')[0])

        alls = set(ids) if alls is None else alls & set(ids)
        print fname, len(ids), len(alls)
    print len(alls)


def std(X):
    xbar = sum(X) / float(len(X))
    return math.sqrt(sum((x - xbar)**2 for x in X)/(len(X) - 1))

def mean(array):
    return float(sum(array))/len(array)