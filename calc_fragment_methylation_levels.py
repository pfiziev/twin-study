import gzip
import multiprocessing
import re
from utils import *


def process((tw, regions)):
#    highpriority()
    print "Twin:", tw
    cRegion = 0
    cChrom = None
    a_regions = {}
    newReg = True
    sitesInNextRegion = []
    for l in gzip.open(datafiles[tw],'r'):
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
                                                  'regEnd'      : regions[cChrom][cRegion][2] })

                    # add sites from the last region that are found in this one, too (due to a possible 4-bases overlap)
                    for sinrMethLevel in sitesInNextRegion:
                        addTo(a_regions[cChrom][-1], 'methLevel', sinrMethLevel)
                        addTo(a_regions[cChrom][-1], 'sites', 1)
                    sitesInNextRegion = []
                    newReg = False


                if cRegion < len(regions[cChrom]) - 1 and regions[cChrom][cRegion+1][1] <= pos <= regions[cChrom][cRegion+1][2]:
                    sitesInNextRegion.append(methLevel)

                addTo(a_regions[cChrom][-1], 'methLevel', methLevel)
                addTo(a_regions[cChrom][-1], 'sites', 1)


    elapsed('annotate regions: '+tw)
    out = open(datafiles[tw] + '.regions', 'w')
    out.writelines(
        "%(chrNo)s\t%(regNo)s\t%(regStart)d\t%(regEnd)d\t%(sites)d\t%(methLevel)f\n"
        % dict(r, **{'chrNo' : chrNo,
                     'methLevel': r['methLevel']/r['sites'] })
                        for chrNo in sorted(regions) for r in a_regions[chrNo]
                            if r['sites'] >= 3 )
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
    multiprocessing.Pool(processes = 4).map(process, [(tw, regions) for tw in datafiles])

