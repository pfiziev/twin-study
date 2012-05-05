import gzip
import json
import os, glob
import cPickle as pickle
import pprint
import re
import datetime
import sys
from utils import *

#twins = dict(   T1 = 'XB002L1',
#                T2 = 'XB002L2',
#                T3 = 'XB002L3',
#                T4 = 'XB002L5',
#                T5 = 'XB002L6',
#                T6 = 'VB015L1',
#                T7 = 'VB015L2',
#                T8 = 'VB015L3',
#                T9 = 'XA004L5',
#                T10 = 'XA004L6' )

if __name__ == '__main__':
#    to_skip = [twins[a] for a in sys.argv[1:]]
#    print 'Skipping', to_skip
#    DATA_DIR = os.path.join(os.getcwd(), 'data')
    sites = None
    data = {}
    stime = datetime.datetime.now()
    for twin_id, fname in datafiles.items():
#        if fname.endswith('.CGmap.gz'): # and not any([skip in fname for skip in to_skip]):

        cdata = {}

        print fname

        current_sites = set()
        for l in gzip.open(fname, 'r'):
            chrNo, nucl, pos, methType, methSubtype, methLevel, methReads, totalReads = filter(None,re.split(r'\s+',l))
            totalReads = int(totalReads)
            if totalReads >= 4 and methType == 'CG':
                key = '%s %s' % (chrNo, pos)
                current_sites.add(key)
                cdata[key] = float(methLevel)
        print "Elapsed:", datetime.datetime.now() - stime
        sites = current_sites if sites is None else sites & current_sites

        data[twin_id] = cdata
        print len(current_sites), len(sites)

    print "STATS:"
#    stats = {"Skipped" : sys.argv[1:],
#             "Total sites" : len(sites) }

    stats = dict((twin_id, sum(data[twin_id][s] for s in sites)/len(sites)) for twin_id in twins)

    for twin_id in twins:
        print "TWIN:%s\tAGE:%d\tMETH_LEVEL:%f\n" % (twin_id, twin_age[twin_id], stats[twin_id])


    logm('Average methylation: \n'+ json.dumps([(twin_id, twin_age[twin_id], stats[twin_id]) for twin_id in twins]))

        
