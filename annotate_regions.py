from itertools import izip
import json
import multiprocessing
from pprint import pformat
from utils import *

def process(tw1, tw2):
    highpriority()

    elapsed("%s %s" % (tw1, tw2))



if __name__ == '__main__':
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
        anno[chrNo].append({'type' : 'gene',
                            'info': geneId,
                            'start' : start,
                            'end' : end
                            })

        exon_starts = map(int, filter(None, fields[9].split(',')))
        exon_ends = map(int, filter(None, fields[10].split(',')))

        for i in xrange(len(exon_starts)):
            anno[chrNo].append({'type' : 'exon',
                                'info':  {'geneId': geneId, 'exon' : i + 1},
                                'start' : exon_starts[i],
                                'end' : exon_ends[i]
                                })
            if i != len(exon_starts) - 1:
                anno[chrNo].append({'type' : 'intron',
                                    'info':  {'geneId': geneId, 'intron' : i + 1},
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
                            'info' : geneId,
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
        anno[chrNo].append({'type'  : 'repeat|%s' % (fields[12] if fields[12] in ['ERV1', 'L1', 'Alu','centr','MaLR','L2','ERVL','ERVK'] else 'Other'),
                            'info'  : '%s|%s|%s'% tuple(fields[10:13]),
                            'start' : int(fields[6]),
                            'end'   : int(fields[7])
                           })
        anno[chrNo].append({'type'  : 'repeat|all',
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

