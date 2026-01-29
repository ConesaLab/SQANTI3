import sys
import csv
import bisect
import itertools
from collections import defaultdict, namedtuple

from src.qc_output import write_isoform_hits
from src.utilities.cupcake.tofu.compare_junctions import compare_junctions

from src.module_logging import qc_logger
from src.qc_classes import myQueryTranscripts
from src.classification_utils import (
    calc_exon_overlap, calc_splicesite_agreement, categorize_full_matches, get_diff_tss_tts,
    calc_overlap, categorize_incomplete_matches, get_gene_diff_tss_tts
)


def transcriptsKnownSpliceSites(isoform_hits_name, refs_1exon_by_chr, refs_exons_by_chr, 
                                start_ends_by_gene, trec, genome_dict, nPolyA):
    """
    This function determines if the isoform hits a known splice site, categorizing it as either
    FSM or ISM.
    :param refs_1exon_by_chr: dict of single exon references (chr -> IntervalTree)
    :param refs_exons_by_chr: dict of multi exon references (chr -> IntervalTree)
    :param trec: id record (genePredRecord) to be compared against reference
    :param genome_dict: dict of genome (chrom --> SeqRecord)
    :param nPolyA: window size to look for polyA
    :return: myQueryTranscripts object that indicates the best reference hit
    """

    # Transcript information for a single query id and comparison with reference.

    # Intra-priming: calculate percentage of "A"s right after the end
    if trec.strand == "+":
        pos_TTS = trec.exonEnds[-1]
        seq_downTTS = str(genome_dict[trec.chrom].seq[pos_TTS:pos_TTS+nPolyA]).upper()
    else: # id on - strand
        pos_TTS = trec.exonStarts[0]
        seq_downTTS = str(genome_dict[trec.chrom].seq[pos_TTS-nPolyA:pos_TTS].reverse_complement()).upper()

    percA = float(seq_downTTS.count('A'))/nPolyA*100

    # Generation of the isoform hit object
    isoform_hit = myQueryTranscripts(isoform = trec.id, chrom = trec.chrom, 
                                     start = trec.txStart, end = trec.txEnd,
                                     strand = trec.strand, length = trec.length, exons = trec.exonCount,
                                     perc_A_downstream_TTS = percA,
                                     seq_A_downstream_TTS = seq_downTTS)


    ##***************************************##
    ########### SPLICED TRANSCRIPTS ###########
    ##***************************************##

    cat_ranking = {'full-splice_match': 5, 'incomplete-splice_match': 4, 'anyKnownJunction': 3, 'anyKnownSpliceSite': 2,
                   'geneOverlap': 1, 'antisense':0.5, '': 0}
    
    # This treats different multi or mono-exonic queries
    if trec.exonCount >= 2:

        hits_by_gene = defaultdict(lambda: [])  # gene --> list of hits
        best_by_gene = {}  # gene --> best isoform_hit

        if trec.chrom in refs_exons_by_chr: # multi-exonic references
            for ref in refs_exons_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
                hits_by_gene[ref.gene].append(ref)
        if trec.chrom in refs_1exon_by_chr: # monoexonic references
            for ref in refs_1exon_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
                hits_by_gene[ref.gene].append(ref)

        if len(hits_by_gene) == 0: return isoform_hit

        for ref_gene, references in hits_by_gene.items(): # Iterate over the matching genes

            # Regenerate the isoform hit object for each gene
            isoform_hit = myQueryTranscripts(isoform = trec.id, chrom = trec.chrom, 
                                            start = trec.txStart, end = trec.txEnd,
                                            strand = trec.strand, length = trec.length, 
                                            exons = trec.exonCount, subcategory= "no_subcategory",
                                            perc_A_downstream_TTS = percA,
                                            seq_A_downstream_TTS = seq_downTTS)

            for ref in references: # Iterate over the matching transcripts for that gene
                exon_overlap = calc_exon_overlap(trec.exons, ref.exons)
                splicesite_hit = calc_splicesite_agreement(trec.exons, ref.exons)   
                
                if trec.strand != ref.strand:
                    # opposite strand, just record it in AS_genes
                    isoform_hit.update({"structural_category": "antisense",
                                        "genes": [ref.gene],
                                        "transcripts": ["novel"],
                                        "q_splicesite_hit": splicesite_hit,
                                        "q_exon_overlap": exon_overlap,
                                        "ref_start": ref.txStart,
                                        "ref_end": ref.txEnd,
                                        "ref_strand":ref.strand,
                                        "ref_exons": ref.exonCount})
                    isoform_hit.AS_genes.add(ref.gene)
                    continue


                
                if ref.exonCount == 1: # mono-exonic reference, handle specially here
                    # TODO: Add a better_than method for the myQueryTranscripts class
                    if exon_overlap > 0 and cat_ranking[isoform_hit.structural_category] < cat_ranking["geneOverlap"]: #CHECK wouldn't it be better with just a number?
                        isoform_hit.update({"structural_category": "geneOverlap",
                                            "subcategory": "mono-exon",
                                            "genes": [ref.gene],
                                            "transcripts": [ref.id],
                                            "ref_length": ref.length,
                                            "ref_exons": ref.exonCount,
                                            "ref_start": ref.txStart,
                                            "ref_end": ref.txEnd,
                                            "ref_strand":ref.strand,
                                            "q_splicesite_hit": splicesite_hit,
                                            "q_exon_overlap": exon_overlap})

                else: # multi-exonic reference
                    match_type = compare_junctions(trec, ref, internal_fuzzy_max_dist=0, max_5_diff=999999, max_3_diff=999999)
                    if match_type not in ('exact', 'subset', 'partial', 'concordant', 'super', 'nomatch'):
                        qc_logger.error(f"Unknown match category {match_type} for isoform {trec.id} against reference {ref.id}!")
                        sys.exit(1)
                    diff_tss, diff_tts = get_diff_tss_tts(trec, ref)

                    # #############################
                    # SQANTI's full-splice_match
                    # #############################

                    if match_type == "exact":
                        subcategory = "multi-exon"
                        if isoform_hits_name:
                            write_isoform_hits(isoform_hits_name,
                                              [trec.id, trec.length, trec.exonCount, ref.id, ref.length, 
                                               ref.exonCount,'FSM', diff_tss, diff_tts])

                        #record hit to isoform hits file
                        # assign as a new hit if
                        # (1) no prev hits yet
                        # (2) this one is better (prev not FSM or is FSM but worse tss/tts)
                        if cat_ranking[isoform_hit.structural_category] < cat_ranking["full-splice_match"] or \
                                                    abs(diff_tss)+abs(diff_tts) < isoform_hit.get_total_diff():
                            subcategory = categorize_full_matches(diff_tss, diff_tts)

                            isoform_hit.update({"structural_category": "full-splice_match",
                                                "subcategory": subcategory,
                                                "genes": [ref.gene],
                                                "transcripts": [ref.id],
                                                "ref_length": ref.length,
                                                "ref_exons": ref.exonCount,
                                                "ref_start": ref.txStart,
                                                "ref_end": ref.txEnd,
                                                "ref_strand": ref.strand,
                                                "diff_to_TSS": diff_tss,
                                                "diff_to_TTS": diff_tts,
                                                "q_splicesite_hit": splicesite_hit,
                                                "q_exon_overlap": exon_overlap})

                    # #######################################################
                    # SQANTI's incomplete-splice_match
                    # (only check if don't already have a FSM match)
                    # #######################################################
                    elif match_type == "subset":
                        subcategory = categorize_incomplete_matches(trec, ref)
                        if isoform_hits_name:
                            write_isoform_hits(isoform_hits_name,
                                              [trec.id, trec.length, trec.exonCount, ref.id, ref.length, 
                                               ref.exonCount, 'ISM', diff_tss, diff_tts])
                        # assign as a new (ISM) hit if
                        # (1) no prev hit
                        # (2) prev hit not as good (is ISM with worse tss/tts or anyKnownSpliceSite)
                        if cat_ranking[isoform_hit.structural_category] < cat_ranking["incomplete-splice_match"] or \
                            (isoform_hit.structural_category=='incomplete-splice_match' and abs(diff_tss)+abs(diff_tts) < isoform_hit.get_total_diff()):
                            isoform_hit.update({"structural_category": "incomplete-splice_match",
                                                "subcategory": subcategory,
                                                "genes": [ref.gene],
                                                "transcripts": [ref.id],
                                                "ref_length": ref.length,
                                                "ref_exons": ref.exonCount,
                                                "ref_start": ref.txStart,
                                                "ref_end": ref.txEnd,
                                                "ref_strand":ref.strand,
                                                "diff_to_TSS": diff_tss,
                                                "diff_to_TTS": diff_tts,
                                                "q_splicesite_hit": splicesite_hit,
                                                "q_exon_overlap": exon_overlap})
                            
                    # #######################################################
                    # Some kind of junction match that isn't ISM/FSM
                    # #######################################################
                    elif match_type in ('partial', 'concordant', 'super'):

                        q_exon_d = abs(trec.exonCount - ref.exonCount)
                    
                        if cat_ranking[isoform_hit.structural_category] < cat_ranking["anyKnownJunction"] or \
                                (isoform_hit.structural_category=='anyKnownJunction' and splicesite_hit > isoform_hit.q_splicesite_hit) or \
                                (isoform_hit.structural_category=='anyKnownJunction' and splicesite_hit == isoform_hit.q_splicesite_hit and exon_overlap > isoform_hit.q_exon_overlap) or \
                                (isoform_hit.structural_category=='anyKnownJunction' and splicesite_hit == isoform_hit.q_splicesite_hit and q_exon_d < abs(trec.exonCount-isoform_hit.ref_exons)):

                            isoform_hit.update({"structural_category": "anyKnownJunction",
                                                "subcategory": "no_subcategory",
                                                "genes": [ref.gene],
                                                "transcripts": ["novel"],
                                                "ref_length": ref.length,
                                                "ref_exons": ref.exonCount,
                                                "ref_start": ref.txStart,
                                                "ref_end": ref.txEnd,
                                                "ref_strand":ref.strand,
                                                "q_splicesite_hit": splicesite_hit,
                                                "q_exon_overlap": exon_overlap})

                    else: # must be nomatch
                        assert match_type == 'nomatch'
                        # at this point, no junction overlap, but may be a single splice site (donor or acceptor) match?
                        # also possibly just exonic (no splice site) overlap
                        if cat_ranking[isoform_hit.structural_category] < cat_ranking["anyKnownSpliceSite"] and splicesite_hit > 0:
                            isoform_hit.update({"structural_category": "anyKnownSpliceSite",
                                                "subcategory": "no_subcategory",
                                                "genes": [ref.gene],
                                                "transcripts": ["novel"],
                                                "ref_length": ref.length,
                                                "ref_exons": ref.exonCount,
                                                "ref_start": ref.txStart,
                                                "ref_end": ref.txEnd,
                                                "ref_strand":ref.strand,
                                                "q_splicesite_hit": splicesite_hit,
                                                "q_exon_overlap": exon_overlap})

                        if isoform_hit.structural_category == "": # still not hit yet, check exonic overlap
                            if cat_ranking[isoform_hit.structural_category] < cat_ranking["geneOverlap"] and exon_overlap > 0:
                                isoform_hit.update({"structural_category": "geneOverlap",
                                                    "subcategory": "no_subcategory",
                                                    "genes": [ref.gene],
                                                    "transcripts": ["novel"],
                                                    "ref_length": ref.length,
                                                    "ref_exons": ref.exonCount,
                                                    "ref_start": ref.txStart,
                                                    "ref_end": ref.txEnd,
                                                    "ref_strand":ref.strand,
                                                    "q_splicesite_hit": splicesite_hit,
                                                    "q_exon_overlap": exon_overlap})

            best_by_gene[ref_gene] = isoform_hit
        # now we have the best transcript hit for each gene, we need to select the best one
        # start with the best scoring one (FSM is best) --> can add other genes if they don't overlap
        geneHitTuple = namedtuple('geneHitTuple', ['score', 'rStart', 'rEnd', 'rGene', 'rStrand', 'iso_hit'])
        best_by_gene = [geneHitTuple(cat_ranking[iso_hit.structural_category], iso_hit.ref_start, iso_hit.ref_end,
                                     ref_gene, iso_hit.ref_strand, iso_hit) for ref_gene, iso_hit in best_by_gene.items()]
        best_by_gene = list(filter(lambda x: x.score > 0, best_by_gene)) # filter out non-matches

        if len(best_by_gene) == 0: # no hit
            return isoform_hit

        # sort matching genes by ranking, allow for multi-gene match as long as they don't overlap
        best_by_gene.sort(key=lambda x: (x.score,
                                         x.iso_hit.q_splicesite_hit+(x.iso_hit.q_exon_overlap)*1./sum(e.end-e.start for e in trec.exons)+calc_overlap(x.rStart,x.rEnd,trec.txStart,trec.txEnd)*1./(x.rEnd-x.rStart)-abs(trec.exonCount-x.iso_hit.ref_exons)), reverse=True)  # sort by (ranking score, overlap)
        isoform_hit = best_by_gene[0].iso_hit
        cur_start, cur_end, cur_strand = best_by_gene[0].rStart, best_by_gene[0].rEnd, best_by_gene[0].rStrand

        for t in best_by_gene[1:]: # Add extra genes if they don't overlap
            if t.score==0: break
            if calc_overlap(cur_start, cur_end, t.rStart, t.rEnd) <= 0 and \
                cur_strand == t.rStrand: # no overlap and in the same strand strands
                isoform_hit.add_gene(t.rGene)
                cur_start, cur_end = min(cur_start, t.rStart), max(cur_end, t.rEnd)

    ##***************************************####
    ########### UNSPLICED TRANSCRIPTS ###########
    ##***************************************####
    else: # single exon isoform
        if trec.chrom in refs_1exon_by_chr:
            for ref in refs_1exon_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
                if ref.strand != trec.strand:
                    #isoform_hit.structural_category = "antisense"
                    # opposite strand, just record it in AS_genes
                    isoform_hit.AS_genes.add(ref.gene)
                    isoform_hit.genes = [ref.gene]
                    continue
                diff_tss, diff_tts = get_diff_tss_tts(trec, ref)

                # see if there's already an existing match AND if so, if this one is better
                if isoform_hit.structural_category in "": # no match so far
                    isoform_hit.update({"structural_category": "full-splice_match",
                                        "subcategory": "mono-exon",
                                        "genes": [ref.gene],
                                        "transcripts": [ref.id],
                                        "ref_length": ref.length,
                                        "ref_exons": ref.exonCount,
                                        "diff_to_TSS": diff_tss,
                                        "diff_to_TTS": diff_tts,
                                        "AS_genes":set()})
                    
                elif abs(diff_tss)+abs(diff_tts) < isoform_hit.get_total_diff():
                    isoform_hit.update({"genes": [ref.gene],
                                        "transcripts": [ref.id],
                                        "ref_length": ref.length,
                                        "ref_exons": ref.exonCount,
                                        "diff_to_TSS": diff_tss,
                                        "diff_to_TTS": diff_tts,
                                        "ref_exons": ref.exonCount,
                                        "AS_genes":set()})
                    
        if isoform_hit.structural_category == "" and trec.chrom in refs_exons_by_chr:
            # no hits to single exon genes, let's see if it hits multi-exon genes
            # (1) if it overlaps with a ref exon and is contained in an exon, we call it ISM
            # (2) else, if it spans one or more introns, we call it NIC by intron retention

            for ref in refs_exons_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
                exon_overlap = calc_exon_overlap(trec.exons, ref.exons)
                if exon_overlap == 0:   # no exonic overlap, skip!
                    continue
                if ref.strand != trec.strand and  exon_overlap > (isoform_hit.q_exon_overlap or 0):
                    # opposite strand, just record it in AS_genes
                    isoform_hit.AS_genes.add(ref.gene)
                    isoform_hit.add_gene(ref.gene)
                    continue

                diff_tss, diff_tts = get_diff_tss_tts(trec, ref)
                flag = False
                
                for e in ref.exons:
                    if e.start <= trec.txStart < trec.txEnd <= e.end:
                        isoform_hit.update({"structural_category": "incomplete-splice_match",
                                            "subcategory": "mono-exon",
                                            "genes": [ref.gene],
                                            "transcripts": [ref.id],
                                            "ref_length": ref.length,
                                            "ref_exons": ref.exonCount,
                                            "diff_to_TSS": diff_tss,
                                            "diff_to_TTS": diff_tts,
                                            "q_exon_overlap": exon_overlap})
                      
                        # this is as good a match as it gets, we can stop the search here
                        get_gene_diff_tss_tts(isoform_hit,trec,start_ends_by_gene)
                        flag = True
                        return isoform_hit
                if flag:
                    continue
                # if we haven't exited here, then ISM hit is not found yet
                # instead check if it's NIC by intron retention
                # but we don't exit here since the next gene could be a ISM hit
                if len(ref.junctions) > 0: # This should always be true since we are checking on multi-exon references
                    for (d,a) in ref.junctions:
                        if trec.txStart < d < a < trec.txEnd:
                            isoform_hit.update({"structural_category": "novel_in_catalog",
                                                "subcategory": "mono-exon_by_intron_retention",
                                                "genes": [ref.gene],
                                                "transcripts": ["novel"],
                                                "ref_length": ref.length,
                                                "ref_exons": ref.exonCount,
                                                "q_exon_overlap": exon_overlap
                                                })
                            get_gene_diff_tss_tts(isoform_hit,trec,start_ends_by_gene)

                            # TODO: Why is this commented out?
                            # return isoform_hit

                # if we get to here, means not ISM, so just add a ref gene and categorize further later
                isoform_hit.add_gene(ref.gene)
                if len(isoform_hit.AS_genes) != 0 and exon_overlap > (isoform_hit.q_exon_overlap or 0):
                    # We found a gene in the correct strand that is longer than the AS gene, so we remove it
                    isoform_hit.AS_genes = set()
                isoform_hit.q_exon_overlap = exon_overlap

            if isoform_hit.structural_category == "novel_in_catalog":
                return isoform_hit
            
    if isoform_hit.genes is not None:
        get_gene_diff_tss_tts(isoform_hit,trec,start_ends_by_gene)
        isoform_hit.genes.sort(key=lambda x: start_ends_by_gene[x]['begin'])

    return isoform_hit


def novelIsoformsKnownGenes(isoform_hit, trec, junctions_by_chr, junctions_by_gene):
    """
    At this point: definitely not FSM or ISM, see if it is NIC, NNC, or fusion
    :return isoform_hit: updated isoforms hit (myQueryTranscripts object)
    """
    def has_intron_retention():
        """
        Checks if there is intron retention in the transcript record (trec).

        This function iterates over the exons of the transcript record and uses binary search to find
        if there are any junctions that indicate intron retention within the exons.

        Returns:
            bool: True if intron retention is detected, False otherwise.
        """
        for e in trec.exons: # TODO:check for strandness
            m = bisect.bisect_left(junctions_by_chr[trec.chrom]['da_pairs'][trec.strand], (e.start, e.end))
            if m < len(junctions_by_chr[trec.chrom]['da_pairs'][trec.strand]) and e.start <= junctions_by_chr[trec.chrom]['da_pairs'][trec.strand][m][0] < junctions_by_chr[trec.chrom]['da_pairs'][trec.strand][m][1] < e.end:
                return True
        return False

    ref_genes = isoform_hit.genes
    # at this point, we have already found matching genes/transcripts
    # hence we do not need to update ref_length or ref_exons
    # or tss_diff and tts_diff (always set to "NA" for non-FSM/ISM matches)
    
    isoform_hit.transcripts = ["novel"]
    if len(ref_genes) == 1:
        # hits exactly one gene, must be either NIC or NNC
        ref_gene_junctions = junctions_by_gene[ref_genes[0]]
        # 1. check if all donors/acceptor sites are known (regardless of which ref gene it came from)
        # 2. check if this query isoform uses a subset of the junctions from the single ref hit
        all_junctions_known = True
        all_junctions_in_hit_ref = True
        for d,a in trec.junctions:
            all_junctions_known = all_junctions_known and (d in junctions_by_chr[trec.chrom]['donors']) and (a in junctions_by_chr[trec.chrom]['acceptors'])
            all_junctions_in_hit_ref = all_junctions_in_hit_ref and ((d,a) in ref_gene_junctions)
        if all_junctions_known:
            isoform_hit.structural_category = "novel_in_catalog"
            if all_junctions_in_hit_ref:
                isoform_hit.subcategory = "combination_of_known_junctions"
            else:
                isoform_hit.subcategory = "combination_of_known_splicesites"
        else:
            isoform_hit.structural_category="novel_not_in_catalog"
            isoform_hit.subcategory = "at_least_one_novel_splicesite"
    elif len(ref_genes) > 1: # see if it is fusion
        # list of a ref junctions from all genes, including potential shared junctions
        # NOTE: some ref genes could be mono-exonic so no junctions
        all_ref_junctions = list(itertools.chain(junctions_by_gene[ref_gene] for ref_gene in ref_genes if ref_gene in junctions_by_gene))

        # (junction index) --> number of refs that have this junction
        junction_ref_hit = dict((i, all_ref_junctions.count(junc)) for i,junc in enumerate(trec.junctions))

        # if the same query junction appears in more than one of the hit references, it is not a fusion
        if max(junction_ref_hit.values()) > 1:
            isoform_hit.structural_category = "moreJunctions"
        else:
            isoform_hit.structural_category = "fusion"
            isoform_hit.subcategory = "mono-exon" if trec.exonCount==1 else "multi-exon"

    if has_intron_retention():
        isoform_hit.subcategory = "intron_retention"

    return isoform_hit

def associationOverlapping(isoform_hit, trec, junctions_by_chr):
    """
    Classify the given isoform based on its overlap with known genes and junctions.

    Parameters:
    isoform_hit (object): An object representing the isoform hit, which will be updated with classification information.
    trec (object): An object representing the transcript record, containing information such as exon count, chromosome, start, and end positions.
    junctions_by_chr (dict): A dictionary where keys are chromosome names and values are dictionaries containing junction information.

    Returns:
    object: The updated isoform_hit object with classification information added.
    """
    # at this point: definitely not FSM or ISM or NIC or NNC
    # possibly (in order of preference assignment):
    #  - antisense  (on opp strand of a known gene)
    #  - genic (overlaps a combination of exons and introns), ignore strand
    #  - genic_intron  (completely within an intron), ignore strand
    #  - intergenic (does not overlap a gene at all), ignore strand
    isoform_hit.update({"structural_category": "intergenic",
                         "transcripts": ["novel"],
                         "subcategory": "mono-exon" if trec.exonCount==1 else "multi-exon"})

    if not isoform_hit.genes:
        # completely no overlap with any genes on the same strand
        # check if it is anti-sense to a known gene, otherwise it's genic_intron or intergenic
        if trec.chrom in junctions_by_chr:
            # no hit even on opp strand
            # see if it is completely contained within a junction
            da_pairs = junctions_by_chr[trec.chrom]['da_tree'].find(trec.txStart, trec.txEnd)
            for junction in da_pairs:
                if junction[0] <= trec.txStart <= trec.txEnd <= junction[1]:
                    isoform_hit.structural_category = "genic_intron"
                    break
        else:
            pass # remain intergenic
    else:
        if len(isoform_hit.AS_genes) == 0:
            isoform_hit.structural_category = "genic"
        else:
            # hits one or more genes on the opposite strand
            isoform_hit.structural_category = "antisense"
            isoform_hit.genes = ["novelGene_{g}_AS".format(g=g) for g in isoform_hit.AS_genes]
    
    return isoform_hit
