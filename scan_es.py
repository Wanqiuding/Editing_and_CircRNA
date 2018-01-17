#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Search exon skipping from junctions.bed generated by tophat.
Output .es file.

Usage:
    python scan_exon_skipping.py -j <junctions.bed> -g <gpe>

Input:
    -j  --junc  junctions.bed generated by tophat.
    -g  --gpe   gene structure annotation in gpe format.

Output format:
exon_skipped    gene    junc_pos    strand  sup_reads   unsup_reads_5
chrX:X-X        XX      chrX:X-X    +/-     int         int
Output format (continue):
unsup_reads_3   exonId_5    exonId_3
int             str         str

Created on Tue Sep 03 14:30:11 2013

@author: yansy
Date: 130904, v1.
Revise: 131002, add exonId_5 and exonId_3 to output.
Revise: 140410, sort juncs after loading so that junctions.bed doesn't have to 
        be pre-sorted. 
Revise: 140413, refind bi-partition search of junctions.

TODO: find exon skipping from bam and gpe.
"""

import sys
import argparse
sys.path.append("/mnt/share/dingwq/bin/Modules_py27")
from genestructure import Bed12, GenePredExt


class Junction(Bed12):
    """
    junction generated by tophat, with 5' donor and 3' acceptor attribute, 
    [donor, acceptor) resembles the intron inside junction, 0-based.
    """
    def __init__(self, line):
        Bed12.__init__(self, line)
        self.donor = self.chromStart + self.blockSizes[0]
        self.acceptor = self.chromStart + self.blockStarts[1]


def find_sub_st(juncLst, pos):
    "find sub for start position of gene structure, 0-based"
    st, ed = 0, len(juncLst)-1
    if juncLst[ed].acceptor < pos:
        return ed+1
    while st < ed:
        md = int((st+ed)/2)
        if juncLst[md].acceptor >= pos and juncLst[md-1].acceptor < pos:
            return md
        elif juncLst[md].acceptor < pos:
            st = md+1
        else:
            ed = md-1
    return ed

def find_sub_ed(juncLst, pos):
    "find sub for end position of gene structure, 0-based"
    st, ed = 0, len(juncLst)-1
    if juncLst[0].donor > pos:
        return 0
    while st < ed:
        md = int((st+ed)/2)
        if juncLst[md].donor <= pos and juncLst[md+1].donor > pos:
            return md
        elif juncLst[md].donor > pos:
            ed = md-1
        else:
            st = md+1
    return ed

def find_juncs(junclist, gpe):
    "find juncs inside gpe, return as a list"
    st_sub = find_sub_st(junclist, gpe.txStart)
    ed_sub = find_sub_ed(junclist, gpe.txEnd)
    return junclist[st_sub:ed_sub+1]

def possible_skipping(gpe):
    """
    Generate all possible exon skipping event. return a list of dict:
    {dn:int, ac:int, strand: +/-, inside_exons:[(st, ed),...],
    exon_skipped:[pos_str], 
    exonId_flanking:[(exonId_5(int), exonId_3(int)), ...]}
    """
    skips = []
    for i, dn in enumerate(gpe.exonEnds):
        for j, ac in enumerate(gpe.exonStarts):
            if j-i > 1:
                exon_skipped = []
                exons = []
                exon_5 = i+1 if gpe.strand=="+" else gpe.exonCount-j
                exon_3 = j+1 if gpe.strand=="+" else gpe.exonCount-i
                for inner in range(i+1, j):
                    exons.append((gpe.exonStarts[inner], gpe.exonEnds[inner]))
                    exon_skipped.append("%s:%d-%d" % (gpe.chrom,
                                                      gpe.exonStarts[inner],
                                                      gpe.exonEnds[inner]))
                skips.append({"dn":dn, "ac":ac, "strand":gpe.strand,
                              "exons":exons,
                              "exon_skipped":exon_skipped,
                              "flank_exons":(exon_5, exon_3)})
    return skips



def search_supporting_junc(skip, juncs):
    "return the junc supporting exon skipping, return None if not found"
    for junc in juncs:
        if junc.donor == skip["dn"] and junc.acceptor == skip["ac"]:
            return junc
    return None

def search_internal_juncs(skip, juncs):
    """
    return a 2-element list of juncs inside the exon skipping dict. 
    first element is the most 5' juncs, 2nd element is the most 3' juncs.
    return [] if not found, 
    if only 1 side available, the other side junc is None.
    """
    junc_5p = None
    junc_3p = None
    for junc in juncs:
        if junc.donor == skip["dn"] and junc.acceptor == skip["exons"][0][0]:
            junc_5p = junc 
        elif junc.acceptor == skip["ac"] and junc.donor == skip["exons"][-1][1]:
            junc_3p = junc
    if skip["strand"] == "-":
        junc_5p, junc_3p = junc_3p, junc_5p
    return [junc_5p, junc_3p] if junc_5p != None or junc_3p != None else None

def skip_event(juncs, gpe):
    "find exon skipping events from juncs in the gpe. Output results to stdout"
    # generate all possible exon skipping event
    exonskips = possible_skipping(gpe)
    # match juncs with exon skipping event
    outLst = []
    for skip in exonskips:
        sup_junc = search_supporting_junc(skip, juncs)
        if not sup_junc:
            continue
        internal_juncs = search_internal_juncs(skip, juncs)
        if not internal_juncs:
            continue
        unsup_reads_5 = internal_juncs[0].score if internal_juncs[0] else 0
        unsup_reads_3 = internal_juncs[1].score if internal_juncs[1] else 0
        jpos=sup_junc.chrom+":"+str(sup_junc.donor)+"-"+str(sup_junc.acceptor)
        outLst.append("\t".join([",".join(skip["exon_skipped"]), 
                                 gpe.name, jpos, gpe.strand,
                                 str(sup_junc.score), str(unsup_reads_5), 
                                 str(unsup_reads_3),
                                 "exon"+str(skip["flank_exons"][0]), 
                                 "exon"+str(skip["flank_exons"][1])]))
    if outLst:
        for line in outLst:
            print line
    else:
        print >> sys.stderr, "No exon skipping junctions for %s" % gpe.name
            

def main(fjunc, fgpe, has_bin):
    juncdict = {}
    for line in fjunc:
        if line.startswith("#"):
            continue
        elif line.startswith("track name="):
            continue
        junc = Junction(line)
        if junc.chrom in juncdict:
            juncdict[junc.chrom].append(junc)
        else:
            juncdict[junc.chrom] = [junc]
    fjunc.close()
    
    # sort juncLst of each chrom by junc position, 1st by donor, 2nd acceptor
    for chrom in juncdict:
        juncdict[chrom].sort(key=lambda junc: junc.acceptor)
        juncdict[chrom].sort(key=lambda junc: junc.donor)
        
    # search each gene's possible exon skipping
    print "\t".join(["#exon_skipped", "gene", "junc_pos", "strand",
                     "sup_reads", "unsup_reads_5", "unsup_reads_3",
                     "exonId_5", "exonId_3"])

    for line in fgpe:
        if line.startswith("#"):
            continue
        gpe = GenePredExt(line, bincolumn=has_bin)
        try:
            juncs = find_juncs(juncdict[gpe.chrom], gpe)
        except KeyError:
            print >> sys.stderr, ("chromosome: %s "
            "not in junctions provided") % gpe.chrom
            continue
        if juncs:
            skip_event(juncs, gpe)
        else:
            print >> sys.stderr, "Found no junction for %s" % gpe.name
    fgpe.close()    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-j","--junc", type=argparse.FileType("r"), 
                        help="junctions.bed generated by tophat.")
    parser.add_argument("-g", "--gpe", type=argparse.FileType("r"), 
                        help="gene structure annotation in gpe format")
    parser.add_argument("-b", "--bin", action='store_true',
                        help="specify if the gpe has bin column")
    args = parser.parse_args()

    main(fjunc=args.junc, fgpe=args.gpe, has_bin=args.bin)