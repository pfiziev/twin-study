import gzip
from itertools import izip
import json
import multiprocessing
from pprint import pformat
import re
from utils import *


if __name__ == '__main__':
    anno = read_annotation()
    


#    anno_out = open(os.path.join(ANNO_DIR,'annotation.info'), 'w')
#    anno_out.write(pformat(anno))
#    anno_out.close()
#    elapsed('save annotation')


    cChrom, cAnno, cPos = None, None, None

    out = open(os.path.join(ANNO_DIR,'common_sites.annotated'), 'w')
    sites = None


    for _, fname in datafiles.items():

        print fname
        c_sites = set()
        for l in gzip.open(fname):
            chrNo, nucl, pos, methType, methSubtype, methLevel, methReads, totalReads = filter(None,re.split(r'\s+',l))
            if int(totalReads) >= 4 and methType == 'CG':
                c_sites.add((int(chrNo), int(pos)))
        sites = c_sites if sites is None else sites & c_sites
        elapsed(fname)

    print "sites in common:", len(sites)
    logm("sites in common: "+ str(len(sites)))
    for chrNo, site_pos in sorted(sites):
        chrNo = 'chr%s' % (str(chrNo) if chrNo <= 22 else ['X','Y','M'][chrNo-23])

        if cChrom != chrNo:
            cChrom = chrNo
            cAnno = anno.get(cChrom, [])
            cPos = 0

        tag = json.dumps([])



        while cPos < len(cAnno) and cAnno[cPos]['end'] < site_pos:
            cPos += 1

        if cPos < len(cAnno):
            iPos = cPos
            regs = []
            while iPos < len(cAnno) and cAnno[iPos]['start'] <= site_pos:
                if site_pos <= cAnno[iPos]['end']:
                    regs.append(cAnno[iPos])
                iPos += 1
            tag = json.dumps(regs)

        out.write("\t".join(map(str, [chrNo, site_pos, tag])) + "\n")




    elapsed('annotating')
    out.close()



#    process([(tw, regions) for tw in datafiles.keys()[:1]][0])
#    multiprocessing.Pool(processes = 5).map(process, [(tw, annotation) for tw in datafiles])

