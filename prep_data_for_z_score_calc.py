import json
import multiprocessing
import re
from utils import *


def process((tw, regions)):
    highpriority()
    print "Twin:", tw
    cRegion = 0
    cChrom = None
    a_regions = {}
    newReg = True
    sitesInNextRegion = []
    for l in open(datafiles[tw],'r').readlines():
        chrNo, nucl, pos, methType, methSubtype, methLevel, methReads, totalReads = filter(None,re.split(r'\s+',l))
        totalReads = int(totalReads)
        methLevel = float(methLevel)
        pos = int(pos)
        if totalReads >= 4 and methType == 'CG':
            if cChrom != chrNo:
                cChrom = chrNo
                cRegion = 0
                a_regions[cChrom] = []
                newReg = True


            while cRegion < len(regions[cChrom]) and regions[cChrom][cRegion][2] < pos:
                cRegion += 1
                newReg = True

            if pos < regions[cChrom][cRegion][1]:
                print "Cannot find region for: ", l
                continue

            if cRegion < len(regions[cChrom]):
                if newReg:
                    a_regions[cChrom].append( {   'regNo'       : regions[cChrom][cRegion][0],
                                                  'regStart'    : regions[cChrom][cRegion][1],
                                                  'regEnd'      : regions[cChrom][cRegion][2],
                                                  'sites'       : []})

                    # add sites from the last region that are found in this one, too (due to a possible 4-bases overlap)
                    for sinrPos, sinrMethLevel in sitesInNextRegion:
                        a_regions[cChrom][-1]['sites'].append({'pos' : sinrPos, 'methLevel' : sinrMethLevel })

                    sitesInNextRegion = []
                    newReg = False


                if cRegion < len(regions[cChrom]) - 1 and regions[cChrom][cRegion+1][1] <= pos <= regions[cChrom][cRegion+1][2]:
                    sitesInNextRegion.append((pos, methLevel))

                a_regions[cChrom][-1]['sites'].append({'pos' : pos, 'methLevel' : methLevel })


    elapsed('annotate regions: '+tw)
    out = open(datafiles[tw] + '.fragments_and_sites', 'w')
    out.writelines(
        "%(chrNo)s\t%(regNo)s\t%(regStart)d\t%(regEnd)d\t%(jsites)s\n"
        % dict(r, **{'chrNo' : chrNo,
                     'jsites': json.dumps(r['sites']) })
                        for chrNo in sorted(regions) for r in a_regions[chrNo]
                            if len(r['sites']) >= 3 )
    out.close()
    elapsed(tw)


if __name__ == '__main__':
    regions = {}
    for l in open(regionsinfo,'r').readlines():
        chrNo, regNo, regStart, regEnd = (lambda (a1,a2,a3,a4): [a1, int(a2), int(a3), int(a4)])(l.split('\t'))
        if chrNo not in regions:
            regions[chrNo] = []
        regions[chrNo].append((regNo, regStart, regEnd))
    elapsed("regions")
#    process([(tw, regions) for tw in datafiles.keys()[:1]][0])
    multiprocessing.Pool(processes = 5).map(process, [(tw, regions) for tw in datafiles])

