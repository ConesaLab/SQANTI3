#!/usr/bin/env python
import os, re, sys, time, subprocess, argparse, pdb
from collections import namedtuple, Counter, defaultdict
from csv import DictReader, DictWriter

from Bio import SeqIO

# Written by Hector del Risco - hdelrisco@ufl.edu

# Definitions
PATSEQLEN = 10    # max size of sequence area to search for match (does not include wiggle)


#### Function to check splice junctions for possible RT switching

#
# Load splice junctions
# Expected file format:
#
# id    junctionNumber  chrom   strand  genomicStartCoord   genomicEndCoord transcriptCoord ...
#
# It returns a hash of all the transcripts with a corresponding array of SpliceJunctions
# (namedtuple is less efficient but it makes the code easier to read - can change to plain array for performance)
# it also returns SJCounts
#
SpliceJunctions = namedtuple("SpliceJunctions", "trans, sjn, chromo, strand, strpos, endpos, transpos, category, startCat, endCat, type")
SJCounts = namedtuple("SJCounts", "trans, sjTotal, sj, knownCanonical, knownNonCanonical, novelCanonical, novelNonCanonical")

FIELDS_RTS = ['isoform', 'junction_number', 'chrom', 'strand', 'genomic_start_coord', 'genomic_end_coord', 'category', 'type', 'exonSeq', 'intronSeq', 'matchLen', 'matchPat', 'mismatch']


def loadSpliceJunctions(filepath):
    """
    Process a splice junction file by SQANTI (see FIELDS_JUNC in sqanti_qc.py)
    :param filepath: the junctions.txt file
    :return: sj_dict (isoform --> junction info), sj_seen_counts ((chr,strand,start,end) --> count of this junction)
    """
    sj_dict = {}
    #sj_type_counts = {'known_canonical':0, 'known_non_canonical':0, 'novel_canonical':0, 'novel_non_canonical':0}
    sj_seen_counts = Counter()

    for rec in DictReader(open(filepath), delimiter='\t'):
        trans = rec['isoform']
        if trans not in sj_dict:
            sj_dict[trans] = []

        # check for unique splice junctions, does not make sense to use duplicates
        sj_pair = (rec['chrom'], rec['strand'], rec['genomic_start_coord'], rec['genomic_end_coord'])
        if sj_pair not in sj_seen_counts:
            sj_seen_counts[sj_pair] = 1
            assert rec['canonical'] in ('canonical', 'non_canonical')
            assert rec['junction_category'] in ('novel', 'known')
            sj_dict[trans].append(SpliceJunctions(trans,
                                                chromo=rec['chrom'],
                                                strand=rec['strand'],
                                                sjn=rec['junction_number'],
                                                strpos=int(rec['genomic_start_coord']),
                                                endpos=int(rec['genomic_end_coord']),
                                                transpos=None,  # ToDo: fix this later
                                                category=rec['junction_category'],
                                                startCat=rec['start_site_category'],
                                                endCat=rec['end_site_category'],
                                                type=rec['canonical']))
            #sj_type_counts[rec['junction_category']+'_'+rec['canonical']] += 1
        else:
            sj_seen_counts[sj_pair] += 1

    return sj_dict, dict(sj_seen_counts)


def checkSJforRTS(sj_dict, genome_dict, wiggle_count, include_category, include_type, min_match, allow_mismatch, output_filename):
    """
    :param sj_dict: dict of (isoform --> junction info)
    :param genome_dict: dict of (chr --> SeqRecord)
    :return: dict of (isoform) -> list of RT junctions. NOTE: dict[isoform] = [] means all junctions are not RT.
    """
    RTS_info_by_isoform = {} # isoform -> list of junction numbers that have RT (ex: 'PB.1.1' --> ['junction_1'])

    wiggle = wiggle_count
    cnt = PATSEQLEN + (2 * wiggle)

    f = open(output_filename, 'w')
    fout = DictWriter(f, fieldnames=FIELDS_RTS, delimiter='\t')
    fout.writeheader()

    for isoform in sj_dict:
        RTS_info_by_isoform[isoform] = []
        #if isoform=='PB.3441.1': pdb.set_trace()
        # process all splice junctions
        for sj in sj_dict[isoform]:
            if (include_type=='c' and sj.type!='canonical') or \
                (include_type=='n' and sj.type!='non_canonical') or \
                (include_category=='n' and sj.category!='novel') or \
                (include_category=='k' and sj.category!='known'):
                continue

            # NOTE: sj.strpos and sj.endpos are both 1-based!!
            # get sequences for pattern and search area
            # the SJ start and end position are positioned at the start/end of the intron
            # ths SJ start is always a lower position than the end regardless of strand
            if sj.strand == "+":
                # we always subtract to get to starting position
                # the count includes 2 wiggles, both ends, so we adjust by 1 wiggle when positioning
                # sequence data on disk: lowpos ----> hipos
                # 5' -----exonSeq(SJstrpos)--------intronSeq(SJendpos) 3'
                _start = sj.strpos - cnt + wiggle - 1
                _end = sj.endpos - cnt + wiggle
                seq_exon = str(genome_dict[sj.chromo].seq[_start:_start+cnt]).upper()
                seq_intron = str(genome_dict[sj.chromo].seq[_end:_end+cnt]).upper()
            else:
                # we are almost on the starting position so just a minor adjustment
                # sequence data on disk: lowpos ----> hipos
                # 3' -----(SJstrpos)intronSeq--------(SJendpos)exonSeq 5'
                _end = sj.strpos - wiggle - 1
                seq_intron = str(genome_dict[sj.chromo][_end:_end+cnt].seq.reverse_complement()).upper()
                _start = sj.endpos - wiggle
                seq_exon = str(genome_dict[sj.chromo][_start:_start+cnt].seq.reverse_complement()).upper()

            # check for RTS repeats and save results to file
            if len(seq_exon) > 0 and len(seq_intron) > 0:
                flag, matchLen, matchPat, mismatch = checkForRepeatPat(seq_exon, seq_intron, min_match, allow_mismatch)
                if flag:
                    RTS_info_by_isoform[isoform].append(sj.sjn)

                    rec = {'isoform': isoform,
                           'junction_number': sj.sjn,
                           'chrom': sj.chromo,
                           'strand': sj.strand,
                           'genomic_start_coord': sj.strpos,
                           'genomic_end_coord': sj.endpos,
                           'category': sj.category,
                           'type': sj.type,
                           'exonSeq': seq_exon,
                           'intronSeq': seq_intron,
                           'matchLen': matchLen,
                           'matchPat': matchPat,
                           'mismatch': mismatch }
                    fout.writerow(rec)

    return RTS_info_by_isoform

# Check for possible RTS
#
# Because at most 1 mismatch is allowed, we can look for exact k-mer matches where k=<min_match>/2
# As soon as we find a matching segment, return True
#
def checkForRepeatPat(seq_exon, seq_intron, min_match, allow_mismatch=True):
    """
    :return: is_RTS (bool), matchLen, matchPattern, mismatch
    """
    #if seq_exon=='TAAATGTACAGG': pdb.set_trace()
    seedsize = int(min_match/2)
    n = len(seq_exon)
    for i in range(n-seedsize+1):
        seed = seq_exon[i:i+seedsize]
        offset = 0
        while True:
            j = seq_intron.find(seed, offset)
            if j >= 0:
                offset = j + 1
                # try to extend the match
                k = seedsize
                while i+k < n and j+k < n and seq_exon[i+k]==seq_intron[j+k]: k += 1
                # right now seq_exon[i:i+k] == seq_intron[j:j+k]
                # we need (min_match - k) more matches on either side
                m = min_match - k
                if (i+k+m <= n) and (j+k+m <= n):
                    flag, mismatch = seq_match(seq_exon[i+k:i+k+m], seq_intron[j+k:j+k+m], allow_mismatch)
                    if flag:
                        return True, k+m, seq_exon[i:i+k+m], mismatch
                if (i-m >= 0) and (j-m >= 0): # try extending mismatch the other way
                    flag, mismatch = seq_match(seq_exon[i-m:i], seq_intron[j-m:j], allow_mismatch)
                    if flag:
                        return True, k+m, seq_exon[i-m:i+k], mismatch
            else:
                break

    return False, None, None, None

#
# Check if sequences match - sequences must have the same length
#
# Note: If mismatch flag is set, will allow 1 mismatch in comparison
#       Regardless of mismacth flag value, indels are not allowed
#
def seq_match(exseq, inseq, allowMismatch):
    """
    Return True if <exseq> and <inseq> are same length and either
    (1) identical OR
    (2) has at most one mismatch (if allowMismatch is True)

    :return: bool, num_mismatch
    """
    if len(exseq)!=len(inseq):
        return False, None
    elif exseq == inseq:
        return True, 0
    elif allowMismatch:
        # allow at most one mismatch
        num_mismatch = 0
        for a,b in zip(exseq, inseq):
            if a!=b:
                if num_mismatch == 1: return False, None  # second mismatch, return False!
                else: num_mismatch += 1
        return True, num_mismatch
    else:
        return False, None


def rts(args, genome_dict):
    #
    # Visualization of where the pattern repeat is searched for using a wiggle (w) of 1
    #
    #       -------------------------                     ---------------------
    #                       w.....cba|w ------- w.....cba|w
    #       -------------------------                     ----------------------
    #                5' exon                intron              3' exon
    #
    # A perfect match is from the end (3') of the first exon, to the end (3') of the first intron
    # This applies to both the + and - strand (- strand is modified by the read function for 5' to '3 and complement)
    # In the case above a perfect match is "abc....." matching "abc....."
    #

    parser = get_parser()
    args = parser.parse_args(args)
    #
    # Check splice junctions for RTS and write results to results file
    #
    absDir = os.path.dirname(os.path.abspath(args.sjFilepath))    
    rts_dir = os.path.join(absDir, "RTS")
    if not os.path.exists(rts_dir):
        os.makedirs(rts_dir)

    rtsResultsFilepath = os.path.join(rts_dir, "sj.rts.results.tsv")

    # load required data
    sjIdx, sjCounts = loadSpliceJunctions(args.sjFilepath)

    # perform RTS analysis
    RTSinfo = checkSJforRTS(sjIdx, genome_dict, args.wiggle_count, args.include_category, args.include_type,
                            args.min_match, args.allow_mismatch, rtsResultsFilepath)

    return RTSinfo

def get_parser():
    # parse command line arguments
    parser = argparse.ArgumentParser(description="Check splice junctions for possible RT switching")
    parser.add_argument('sjFilepath', type=str, help='file with splice junction information')
    parser.add_argument('mmfaFilepath', type=str, help='path to reference genome')
    parser.add_argument("-m", "--min_match", type=int, default=8, choices=list(range(4, 11)), help="Minimum number of bases required to match. Default: 8")
    parser.add_argument("-a", "--allow_mismatch", default=False, action="store_true", help="Specify to allow 1 base mismatch in sequences (indels are not allowed)")
    parser.add_argument("-w", "--wiggle_count", type=int, default=1, choices=list(range(0, 4)), help="Number of bases allowed to wiggle on each side of ideal RTS sequence location. Default: 1")
    parser.add_argument("-t", "--include_type", default='a', choices=['a', 'c', 'n'], help="Type of splice junctions to include (a for all, c for canonical, and n for non-canonical). Default: a")
    parser.add_argument("-c", "--include_category", default='a', choices=['a', 'n', 'k'], help="Category of splice junctions to include (a for all, n for novel, and k for known). Default: a")
    parser.add_argument("-v", "--version", help="Display program version number", action='version', version='%(prog)s 0.1')

    return parser

if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()

    print("Reading genome fasta into dict...", file=sys.stderr)
    genome_dict = dict((r.id, r) for r in SeqIO.parse(open(args.mmfaFilepath), 'fasta'))
