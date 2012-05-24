from itertools import izip
import json
import multiprocessing
from pprint import pformat
from utils import *

def process(tw1, tw2):
    highpriority()

    elapsed("%s %s" % (tw1, tw2))



if __name__ == '__main__':
    anno = read_annotation()

#    anno_out = open(os.path.join(ANNO_DIR,'annotation.info'), 'w')
#    anno_out.write(pformat(anno))
#    anno_out.close()
#    elapsed('save annotation')


    cChrom, cAnno, cPos = None, None, None

    out = open(os.path.join(ANNO_DIR,'RRBS_mapable_regions.info.annotated'), 'w')
    for l in open(os.path.join(ANNO_DIR,'RRBS_mapable_regions.info'), 'r'):
        chrNo, regId, regStart, regEnd = map(int, l.split('\t'))
        chrNo = 'chr%s' % (str(chrNo) if chrNo <= 22 else ['X','Y','M'][chrNo-23])

        if cChrom != chrNo:
            cChrom = chrNo
            cAnno = anno.get(cChrom, [])
            cPos = 0

        tag = json.dumps([])

        def total_overlap(s1, e1, s2, e2):
            return s1 <= s2 and e1 >= e2 or s2 <= s1 and e2 >= e1

        def partial_overlap(s1, e1, s2, e2):
            return s1 <= e2 and s2 <= e1 and min(e1, e2) - max(s1, s2) >= 10


        while cPos < len(cAnno) and cAnno[cPos]['end'] < regStart:
            cPos += 1

        if cPos < len(cAnno):
            iPos = cPos
            regs = []
            while iPos < len(cAnno) and cAnno[iPos]['start'] <= regEnd:
                if partial_overlap(cAnno[iPos]['start'], cAnno[iPos]['end'], regStart, regEnd):
                    regs.append(cAnno[iPos])
                iPos += 1
            tag = json.dumps(regs)
            
        out.write("\t".join(map(str, [chrNo, regId, regStart, regEnd, tag])) + "\n")




    elapsed('annotating')
    out.close()

        

#    process([(tw, regions) for tw in datafiles.keys()[:1]][0])
#    multiprocessing.Pool(processes = 5).map(process, [(tw, annotation) for tw in datafiles])

