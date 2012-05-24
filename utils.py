from itertools import product, chain
import os
import datetime
import math
import re


twinpairs = {}
GAY_TWINS = []
HETERO_TWINS = []
twin_age = {}
_twins = {}

_path = os.path.split(__file__)[0]
pwd = lambda p: os.path.join(_path, p)


twin_info = open(pwd('CGmap/twin_info'))
for line in twin_info:
    if line.startswith('#') or line.strip() == '': continue
    _twin_pair_id, _twin_id, _gay, _filename, _age = line.strip().split('\t')

    twin_age[_twin_id] = int(_age)
    (GAY_TWINS if _gay == '1' else HETERO_TWINS).append(_twin_id)
    if _twin_pair_id not in twinpairs: twinpairs[_twin_pair_id] = []
    twinpairs[_twin_pair_id].append(_twin_id)
    _twins[_twin_id] = _filename

twin_info.close()
twinpairs = map(tuple, twinpairs.values())

#_twins = dict(  T1 = 'XB002L1',
#                T2 = 'XB002L2',
#                T3 = 'XB002L3',
#                T4 = 'XB002L5',
#                T5 = 'XB002L6',
#                T6 = 'VB015L1',
#                T7 = 'VB015L2',
#                T8 = 'VB015L3',
#                T9 = 'XA004L5',
#                T10 = 'XA004L6' )
#
#
#twinpairs = [('T7','T8'),
#             ('T5','T6'),
#             ('T1','T2'),
#             ('T9','T10'),
#             ('T3','T4')]
#
## True if gay, False if hetero
#GAY_TWINS    = frozenset(['T1', 'T3', 'T5', 'T7', 'T9' ])
#HETERO_TWINS = frozenset(['T2', 'T4', 'T6', 'T8', 'T10'])


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
#
#twin_age = {'T7'  : 21,
#            'T8'  : 21,
#            'T5'  : 30,
#            'T6'  : 30,
#            'T1'  : 38,
#            'T2'  : 38,
#            'T9'  : 45,
#            'T10' : 45,
#            'T3'  : 55,
#            'T4'  : 55}





#twin_random = {'T1': 23,
# 'T10': 33,
# 'T2': 37,
# 'T3': 31,
# 'T4': 55,
# 'T5': 38,
# 'T6': 36,
# 'T7': 29,
# 'T8': 49,
# 'T9': 47}



twins = sorted(twin_age, key = lambda k: twin_age[k])


OUTPUT_DIR = pwd('output')
def outd(filename):
    return os.path.join(OUTPUT_DIR, filename)


DATA_DIR = pwd('CGmap')

_datafiles = [os.path.join(DATA_DIR, fname) for fname in os.listdir(DATA_DIR) if fname.endswith('.CGmap.gz')]
datafiles = dict((key, filter(lambda x: _twins[key] in x, _datafiles)[0]) for key in _twins)
file2twin = dict((fname, twin_id) for twin_id, fname in datafiles.items())



stime = datetime.datetime.now()
def elapsed(msg = None):
    print "[%s]" % msg if msg is not None else "+", "Elapsed: " , datetime.datetime.now() - stime


ANNO_DIR = pwd('annotation')

regionsinfo = os.path.join(ANNO_DIR, "RRBS_mapable_regions.info")

#OUTPUT_DIR = os.path.join(os.getcwd(),'output')


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


logf = open(os.path.join(DATA_DIR, 'log'), 'a')
def logm(message):
    logf.write('[%s] %s\n' % (str(datetime.datetime.now()), message))



def read_annotation():
    anno = {}
    print "reading refGene"
    for l in open(os.path.join(ANNO_DIR, 'refGene.txt'), 'r'):
        fields = l.split('\t')
        chrNo = fields[2]
        geneId = fields[1]
        if chrNo not in anno:
            anno[chrNo] = []
        start = int(fields[4])
        end = int(fields[5])
        geneSymbol = fields[12]
        anno[chrNo].append({'type' : 'gene',
                            'info': {'geneId': geneId, 'geneSymbol' : geneSymbol},
                            'geneSymbol' : geneSymbol,
                            'start' : start,
                            'end' : end
        })

        exon_starts = map(int, filter(None, fields[9].split(',')))
        exon_ends = map(int, filter(None, fields[10].split(',')))

        for i in xrange(len(exon_starts)):
            anno[chrNo].append({'type' : 'exon',
                                'info':  {'geneId': geneId, 'exon' : i + 1, 'geneSymbol' : geneSymbol},
                                'start' : exon_starts[i],
                                'end' : exon_ends[i]
            })
            if i != len(exon_starts) - 1:
                anno[chrNo].append({'type' : 'intron',
                                    'info':  {'geneId': geneId, 'intron' : i + 1, 'geneSymbol' : geneSymbol},
                                    'start' : exon_ends[i] + 1,
                                    'end' : exon_starts[i + 1] - 1
                })





        if fields[3] == '+':
            promotor_start = start - 3700
            promotor_end = promotor_start + 4000
        else:
            # for minus strand keep start < end (reverse start and end)
            promotor_end = end + 3700
            promotor_start = promotor_end - 4000

        anno[chrNo].append({'type' : 'promotor',
                            'info' : {'geneId': geneId, 'geneSymbol' : geneSymbol},
                            'start' : promotor_start,
                            'end' : promotor_end
        })



    print "reading repeats info"
    for l in open(os.path.join(ANNO_DIR, 'repeats.txt'), 'r'):
        if l.startswith('#bin'): continue
        fields = l.split('\t')
        chrNo = fields[5]
        if chrNo not in anno:
            anno[chrNo] = []
        anno[chrNo].append({'type'  : 'repeat|%s' % (fields[12] if fields[12] in ['ERV1', 'L1', 'Alu','centr','MaLR','L2',
                                                                                  'ERVL','ERVK','Simple_repeat','MIR','MER1_type',
                                                                                  'Low_complexity','L2'] else 'Other'),
                            'info'  : '%s|%s|%s'% tuple(fields[10:13]),
                            'start' : int(fields[6]),
                            'end'   : int(fields[7])
        })
        anno[chrNo].append({'type'  : 'repeat|All',
                            'info'  : '%s|%s|%s'% tuple(fields[10:13]),
                            'start' : int(fields[6]),
                            'end'   : int(fields[7])
        })

    print "reading CpG island info"
    for l in open(os.path.join(ANNO_DIR, 'cpgIslandExt.txt'), 'r'):
        fields = l.split('\t')
        chrNo = fields[0]
        if chrNo not in anno:
            anno[chrNo] = []
        anno[chrNo].append({'type'  : 'CpG island',
                            'start' : int(fields[1]),
                            'end'   : int(fields[2])
        })


    elapsed("reading annotation")
    for chrNo in anno:
        anno[chrNo] = sorted(anno[chrNo], key = lambda x: x['start'])
    elapsed('sorting')


    print sum(len(anno[c]) for c in anno)
    return anno