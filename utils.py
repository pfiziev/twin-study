from itertools import product, chain
import os
import datetime


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

twins = _twins.keys()

twinpairs = [('T7','T8'),
             ('T5','T6'),
             ('T1','T2'),
             ('T9','T10'),
             ('T3','T4')]

random_twin_pairs = [pair for pair in product([tp[0] for tp in twinpairs],[tp[1] for tp in twinpairs]) if pair not in twinpairs]

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


def overlap(*args):
    alls = None
    for fname in args:
        ids = [l.split('\t')[0] for l in open(fname,'r')]

        alls = set(ids) if alls is None else alls & set(ids)
        print fname, len(ids), len(alls)
    print len(alls)