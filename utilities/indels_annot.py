#!/usr/bin/env python

import pysam
from collections import defaultdict, Counter
from csv import DictReader, DictWriter
from bx.intervals import Interval

"""
pysam cigar type:
http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples

M	BAM_CMATCH	0
I	BAM_CINS	1
D	BAM_CDEL	2
N	BAM_CREF_SKIP	3
S	BAM_CSOFT_CLIP	4
H	BAM_CHARD_CLIP	5
P	BAM_CPAD	6
=	BAM_CEQUAL	7
X	BAM_CDIFF	8
B	BAM_CBACK	9
"""

MAX_DIST_FROM_JUNC = 10 # maximum distance from junction in which we label "indel near junction"
CIGAR_TYPE_LIST = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B']
FIELDS_INDEL = ['isoform', 'indelStart', 'indelEnd', 'nt', 'nearJunction', 'junctionStart', 'junctionEnd', 'indelType']

def calc_indels_from_sam(samFile):
    """
    Given an aligned SAM file, calculate indel statistics.
    :param samFile: aligned SAM file
    :return: indelsJunc (dict of pbid --> list of junctions near indel), indelsTotal (dict of pbid --> total indels count)
    """
    sam = pysam.AlignmentFile(samFile, "r")
    out_file = samFile[:samFile.rfind('.')]+"_indels.txt"
    fhandle = open(out_file, "w")
    fout = DictWriter(fhandle, fieldnames=FIELDS_INDEL, delimiter='\t')
    fout.writeheader()

    indelsJunc = defaultdict(lambda: [])
    indelsTotal = Counter()


    for read in sam.fetch():
        if read.is_unmapped:
            continue
        cigarLine = read.cigar
        ## reading splice junctions and storing information
        pos_start = read.pos # 0-based start
        spliceSites = []  # list of splice junctions (Interval(donor, acceptor))

        for (cigarType,cigarLength) in cigarLine:
            if CIGAR_TYPE_LIST[cigarType] in ('M', 'D', 'N', 'P', 'B'):
                pos_end = pos_start + cigarLength # 1-based end
                if (CIGAR_TYPE_LIST[cigarType] == 'N'): # skip (intron)
                    spliceSites.append(Interval(pos_start, pos_end))
                pos_start = pos_end

        ## reading indels, comparing with splice junctions and writing information
        pos_start = read.pos # 0-based start

        for (cigarType,cigarLength) in cigarLine:
            if CIGAR_TYPE_LIST[cigarType] in ('M', 'D', 'N', 'P', 'B'):
                pos_end = pos_start + cigarLength # 1-based end

            if CIGAR_TYPE_LIST[cigarType] in ("I", "D"): # insertion or deletion
                pos_indel = pos_start  # 0-based
                pos_end_indel = pos_start+1 if CIGAR_TYPE_LIST[cigarType]=='I' else pos_end # 1-based
                spliceSitesNearIndel = []
                name = str(read.query_name).split("|")[0]

                # indels in the sequence
                indelsTotal[name] += 1

                # indels near spliceSties
                for sj in spliceSites:
                    if abs(pos_indel-sj.start) < MAX_DIST_FROM_JUNC or abs(pos_indel-sj.end+1) < MAX_DIST_FROM_JUNC or \
                       abs(pos_end_indel-1-sj.start) < MAX_DIST_FROM_JUNC or abs(pos_end_indel-sj.end) < MAX_DIST_FROM_JUNC:
                        spliceSitesNearIndel.append(sj)

                rec = {'isoform': name,
                       'indelStart': pos_indel + 1, # make start 1-based
                       'indelEnd': pos_end_indel,
                       'nt': cigarLength,
                       'nearJunction': "FALSE",
                       'junctionStart': 'NA',
                       'junctionEnd': 'NA',
                       'indelType': 'insertion' if CIGAR_TYPE_LIST[cigarType]=='I' else 'deletion'}
                if len(spliceSitesNearIndel)==0:
                    fout.writerow(rec)
                else:
                    rec['nearJunction'] = 'TRUE'
                    for sj in spliceSitesNearIndel:
                        rec['junctionStart'] = sj.start + 1  # make start now 1-based
                        rec['junctionEnd'] = sj.end          # end is already 1-based
                        fout.writerow(rec)
                        indelsJunc[name].append(sj)

            if CIGAR_TYPE_LIST[cigarType] in ('M', 'D', 'N', 'P', 'B'):
                pos_start = pos_end

    sam.close()
    fhandle.close()
    return dict(indelsJunc), indelsTotal


if __name__ == "__main__":
    import sys
    calc_indels_from_sam(sys.argv[1])