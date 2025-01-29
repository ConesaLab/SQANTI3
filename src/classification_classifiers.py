import csv
import bisect
import itertools
from collections import defaultdict, namedtuple

from .utilities.cupcake.tofu.compare_junctions import compare_junctions

from .qc_classes import myQueryTranscripts
from .classification_utils import (
    calc_exon_overlap, calc_splicesite_agreement, get_diff_tss_tts,
    calc_overlap, categorize_incomplete_matches, get_gene_diff_tss_tts
)


def transcriptsKnownSpliceSites(isoform_hits_name, refs_1exon_by_chr, refs_exons_by_chr, start_ends_by_gene, trec, genome_dict, nPolyA):
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
    isoform_hit = myQueryTranscripts(id=trec.id, tts_diff="NA", tss_diff="NA",\
                                    num_exons=trec.exonCount,
                                    length=trec.length,
                                    str_class="", \
                                    chrom=trec.chrom,
                                    strand=trec.strand, \
                                    subtype="no_subcategory",\
                                    percAdownTTS=str(percA),\
                                    seqAdownTTS=seq_downTTS)

    ##***************************************##
    ########### SPLICED TRANSCRIPTS ###########
    ##***************************************##

    cat_ranking = {'full-splice_match': 5, 'incomplete-splice_match': 4, 'anyKnownJunction': 3, 'anyKnownSpliceSite': 2,
                   'geneOverlap': 1, '': 0}
    
    # This treats different multi or mono-exonic queries
    if trec.exonCount >= 2:

        hits_by_gene = defaultdict(lambda: [])  # gene --> list of hits
        best_by_gene = {}  # gene --> best isoform_hit

        if trec.chrom in refs_exons_by_chr:
            for ref in refs_exons_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
                hits_by_gene[ref.gene].append(ref)
        if trec.chrom in refs_1exon_by_chr:
            for ref in refs_1exon_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
                hits_by_gene[ref.gene].append(ref)

        if len(hits_by_gene) == 0: return isoform_hit

        for ref_gene in sorted(hits_by_gene):
            isoform_hit = myQueryTranscripts(id=trec.id, tts_diff="NA", tss_diff="NA", \
                                             num_exons=trec.exonCount,
                                             length=trec.length,
                                             str_class="", \
                                             chrom=trec.chrom,
                                             strand=trec.strand, \
                                             subtype="no_subcategory", \
                                             percAdownTTS=str(percA), \
                                             seqAdownTTS=seq_downTTS)

            for ref in hits_by_gene[ref_gene]:
                if trec.strand != ref.strand:
                    # opposite strand, just record it in AS_genes
                    isoform_hit.AS_genes.add(ref.gene)
                    continue

                if ref.exonCount == 1: # mono-exonic reference, handle specially here
                    # TODO: calc exon overlap is called more than twice. Create a variable to store the value 
                    if calc_exon_overlap(trec.exons, ref.exons) > 0 and cat_ranking[isoform_hit.str_class] < cat_ranking["geneOverlap"]: #CHECK wouldn't it be better with just a number?
                        isoform_hit = myQueryTranscripts(trec.id, "NA", "NA", trec.exonCount, trec.length,
                                                            "geneOverlap",
                                                             subtype="mono-exon",
                                                             chrom=trec.chrom,
                                                             strand=trec.strand,
                                                             genes=[ref.gene],
                                                             transcripts=[ref.id],
                                                             refLen=ref.length,
                                                             refExons=ref.exonCount,
                                                             refStart=ref.txStart,
                                                             refEnd=ref.txEnd,
                                                             q_splicesite_hit=0,
                                                             q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons),
                                                             percAdownTTS=str(percA),
                                                             seqAdownTTS=seq_downTTS)

                else: # multi-exonic reference
                    match_type = compare_junctions(trec, ref, internal_fuzzy_max_dist=0, max_5_diff=999999, max_3_diff=999999)
                    #TODO: Error handling in logs
                    if match_type not in ('exact', 'subset', 'partial', 'concordant', 'super', 'nomatch'):
                        raise Exception("Unknown match category {0}!".format(match_type))

                    diff_tss, diff_tts = get_diff_tss_tts(trec, ref)
                    #has_overlap = gene_overlap(isoform_hit.genes[-1], ref.gene) if len(isoform_hit.genes) >= 1 else Fals
                    # #############################
                    # SQANTI's full-splice_match
                    # #############################
                    if match_type == "exact":
                        subtype = "multi-exon"
                        if isoform_hits_name:
                            with open(isoform_hits_name+'_tmp', 'a') as out_file:
                                tsv_writer = csv.writer(out_file, delimiter='\t')
                                tsv_writer.writerow([trec.id, trec.length, trec.exonCount, ref.id, ref.length, ref.exonCount,
                                                      'FSM', diff_tss, diff_tts])
                        #record hit to isoform hits file
                        # assign as a new hit if
                        # (1) no prev hits yet
                        # (2) this one is better (prev not FSM or is FSM but worse tss/tts)
                        if cat_ranking[isoform_hit.str_class] < cat_ranking["full-splice_match"] or \
                                                    abs(diff_tss)+abs(diff_tts) < isoform_hit.get_total_diff():
                            # subcategory for matching 5' and matching 3'
                            if abs(diff_tss) <= 50 and abs(diff_tts) <= 50:
                                    subtype = 'reference_match'
                            # subcategory for matching 5' and non-matching 3'
                            if abs(diff_tss) <= 50 and abs(diff_tts) > 50:
                                subtype = 'alternative_3end'
                            # subcategory for matching 3' and non-matching 5'
                            if abs(diff_tss) > 50 and abs(diff_tts) <= 50:
                                subtype = 'alternative_5end'
                            # subcategory for non-matching 3' and non-matching 5'
                            if abs(diff_tss) > 50 and abs(diff_tts) > 50:
                                subtype = 'alternative_3end5end'
                            isoform_hit = myQueryTranscripts(trec.id, diff_tss, diff_tts, trec.exonCount, trec.length,
                                                              str_class="full-splice_match",
                                                              subtype=subtype,
                                                              chrom=trec.chrom,
                                                              strand=trec.strand,
                                                              genes=[ref.gene],
                                                              transcripts=[ref.id],
                                                              refLen = ref.length,
                                                              refExons= ref.exonCount,
                                                              refStart=ref.txStart,
                                                              refEnd=ref.txEnd,
                                                              q_splicesite_hit=calc_splicesite_agreement(trec.exons, ref.exons),
                                                              q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons),
                                                              percAdownTTS=str(percA),
                                                              seqAdownTTS=seq_downTTS)
                    # #######################################################
                    # SQANTI's incomplete-splice_match
                    # (only check if don't already have a FSM match)
                    # #######################################################
                    elif match_type == "subset":
                        subtype = categorize_incomplete_matches(trec, ref)
                        if isoform_hits_name:
                            with open(isoform_hits_name+'_tmp', 'a') as out_file:
                                tsv_writer = csv.writer(out_file, delimiter='\t')
                                tsv_writer.writerow([trec.id, trec.length, trec.exonCount, ref.id, ref.length,
                                                     ref.exonCount, 'ISM', diff_tss, diff_tts])
                        # assign as a new (ISM) hit if
                        # (1) no prev hit
                        # (2) prev hit not as good (is ISM with worse tss/tts or anyKnownSpliceSite)
                        if cat_ranking[isoform_hit.str_class] < cat_ranking["incomplete-splice_match"] or \
                            (isoform_hit.str_class=='incomplete-splice_match' and abs(diff_tss)+abs(diff_tts) < isoform_hit.get_total_diff()):
                            isoform_hit = myQueryTranscripts(trec.id, diff_tss, diff_tts, trec.exonCount, trec.length,
                                                             str_class="incomplete-splice_match",
                                                             subtype=subtype,
                                                             chrom=trec.chrom,
                                                             strand=trec.strand,
                                                             genes=[ref.gene],
                                                             transcripts=[ref.id],
                                                             refLen = ref.length,
                                                             refExons= ref.exonCount,
                                                             refStart=ref.txStart,
                                                             refEnd=ref.txEnd,
                                                             q_splicesite_hit=calc_splicesite_agreement(trec.exons, ref.exons),
                                                             q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons),
                                                             percAdownTTS=str(percA),
                                                             seqAdownTTS=seq_downTTS)
                    # #######################################################
                    # Some kind of junction match that isn't ISM/FSM
                    # #######################################################
                    elif match_type in ('partial', 'concordant', 'super'):
                        q_sp_hit = calc_splicesite_agreement(trec.exons, ref.exons)
                        q_ex_overlap = calc_exon_overlap(trec.exons, ref.exons)
                        q_exon_d = abs(trec.exonCount - ref.exonCount)
                        if cat_ranking[isoform_hit.str_class] < cat_ranking["anyKnownJunction"] or \
                                (isoform_hit.str_class=='anyKnownJunction' and q_sp_hit > isoform_hit.q_splicesite_hit) or \
                                (isoform_hit.str_class=='anyKnownJunction' and q_sp_hit==isoform_hit.q_splicesite_hit and q_ex_overlap > isoform_hit.q_exon_overlap) or \
                                (isoform_hit.str_class=='anyKnownJunction' and q_sp_hit==isoform_hit.q_splicesite_hit and q_exon_d < abs(trec.exonCount-isoform_hit.refExons)):
                            isoform_hit = myQueryTranscripts(trec.id, "NA", "NA", trec.exonCount, trec.length,
                                                             str_class="anyKnownJunction",
                                                             subtype="no_subcategory",
                                                             chrom=trec.chrom,
                                                             strand=trec.strand,
                                                             genes=[ref.gene],
                                                             transcripts=["novel"],
                                                             refLen=ref.length,
                                                             refExons=ref.exonCount,
                                                             refStart=ref.txStart,
                                                             refEnd=ref.txEnd,
                                                             q_splicesite_hit=calc_splicesite_agreement(trec.exons, ref.exons),
                                                             q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons),
                                                             percAdownTTS=str(percA),
                                                             seqAdownTTS=seq_downTTS)
                    else: # must be nomatch
                        assert match_type == 'nomatch'
                        # at this point, no junction overlap, but may be a single splice site (donor or acceptor) match?
                        # also possibly just exonic (no splice site) overlap
                        if cat_ranking[isoform_hit.str_class] < cat_ranking["anyKnownSpliceSite"] and calc_splicesite_agreement(trec.exons, ref.exons) > 0:
                            isoform_hit = myQueryTranscripts(trec.id, "NA", "NA", trec.exonCount, trec.length,
                                                             str_class="anyKnownSpliceSite",
                                                             subtype="no_subcategory",
                                                             chrom=trec.chrom,
                                                             strand=trec.strand,
                                                             genes=[ref.gene],
                                                             transcripts=["novel"],
                                                             refLen=ref.length,
                                                             refExons=ref.exonCount,
                                                             refStart=ref.txStart,
                                                             refEnd=ref.txEnd,
                                                             q_splicesite_hit=calc_splicesite_agreement(trec.exons, ref.exons),
                                                             q_exon_overlap=calc_exon_overlap(trec.exons,
                                                                                              ref.exons),
                                                             percAdownTTS=str(percA),
                                                             seqAdownTTS=seq_downTTS)

                        if isoform_hit.str_class=="": # still not hit yet, check exonic overlap
                            if cat_ranking[isoform_hit.str_class] < cat_ranking["geneOverlap"] and calc_exon_overlap(trec.exons, ref.exons) > 0:
                                isoform_hit = myQueryTranscripts(trec.id, "NA", "NA", trec.exonCount, trec.length,
                                                                 str_class="geneOverlap",
                                                                 subtype="no_subcategory",
                                                                 chrom=trec.chrom,
                                                                 strand=trec.strand,
                                                                 genes=[ref.gene],
                                                                 transcripts=["novel"],
                                                                 refLen=ref.length,
                                                                 refExons=ref.exonCount,
                                                                 refStart=ref.txStart,
                                                                 refEnd=ref.txEnd,
                                                                 q_splicesite_hit=calc_splicesite_agreement(trec.exons, ref.exons),
                                                                 q_exon_overlap=calc_exon_overlap(trec.exons, ref.exons),
                                                                 percAdownTTS=str(percA),
                                                                 seqAdownTTS=seq_downTTS)

            best_by_gene[ref_gene] = isoform_hit
        # now we have best_by_gene:
        # start with the best scoring one (FSM is best) --> can add other genes if they don't overlap
  
        geneHitTuple = namedtuple('geneHitTuple', ['score', 'rStart', 'rEnd', 'rGene', 'iso_hit'])
        best_by_gene = [geneHitTuple(cat_ranking[iso_hit.str_class],iso_hit.refStart,iso_hit.refEnd,ref_gene,iso_hit) for ref_gene,iso_hit in best_by_gene.items()]
        best_by_gene = list(filter(lambda x: x.score > 0, best_by_gene))
        if len(best_by_gene) == 0: # no hit
            return isoform_hit

        # sort matching genes by ranking, allow for multi-gene match as long as they don't overlap
        best_by_gene.sort(key=lambda x: (x.score,x.iso_hit.q_splicesite_hit+(x.iso_hit.q_exon_overlap)*1./sum(e.end-e.start for e in trec.exons)+calc_overlap(x.rStart,x.rEnd,trec.txStart,trec.txEnd)*1./(x.rEnd-x.rStart)-abs(trec.exonCount-x.iso_hit.refExons)), reverse=True)  # sort by (ranking score, overlap)
        isoform_hit = best_by_gene[0].iso_hit
        cur_start, cur_end = best_by_gene[0].rStart, best_by_gene[0].rEnd
        for t in best_by_gene[1:]:
            if t.score==0: break
            if calc_overlap(cur_start, cur_end, t.rStart, t.rEnd) <= 0:
                isoform_hit.genes.append(t.rGene)
                cur_start, cur_end = min(cur_start, t.rStart), max(cur_end, t.rEnd)

    ##***************************************####
    ########### UNSPLICED TRANSCRIPTS ###########
    ##***************************************####
    else: # single exon id
        if trec.chrom in refs_1exon_by_chr:
            for ref in refs_1exon_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
                if ref.strand != trec.strand:
                    # opposite strand, just record it in AS_genes
                    isoform_hit.AS_genes.add(ref.gene)
                    continue
                diff_tss, diff_tts = get_diff_tss_tts(trec, ref)

                # see if there's already an existing match AND if so, if this one is better
                if isoform_hit.str_class == "": # no match so far
                    isoform_hit = myQueryTranscripts(trec.id, diff_tss, diff_tts, trec.exonCount, trec.length, "full-splice_match",
                                                            subtype="mono-exon",
                                                            chrom=trec.chrom,
                                                            strand=trec.strand,
                                                            genes=[ref.gene],
                                                            transcripts=[ref.id],
                                                            refLen=ref.length,
                                                            refExons = ref.exonCount,
                                                            percAdownTTS=str(percA),
                                                            seqAdownTTS=seq_downTTS)
                elif abs(diff_tss)+abs(diff_tts) < isoform_hit.get_total_diff():
                    isoform_hit.modify(ref.id, ref.gene, diff_tss, diff_tts, ref.length, ref.exonCount)

        if isoform_hit.str_class == "" and trec.chrom in refs_exons_by_chr:
            # no hits to single exon genes, let's see if it hits multi-exon genes
            # (1) if it overlaps with a ref exon and is contained in an exon, we call it ISM
            # (2) else, if it spans one or more introns, we call it NIC by intron retention
            for ref in refs_exons_by_chr[trec.chrom].find(trec.txStart, trec.txEnd):
                if calc_exon_overlap(trec.exons, ref.exons) == 0:   # no exonic overlap, skip!
                    continue
                if ref.strand != trec.strand:
                    # opposite strand, just record it in AS_genes
                    isoform_hit.AS_genes.add(ref.gene)
                    continue
                diff_tss, diff_tts = get_diff_tss_tts(trec, ref)
                flag = False
                for e in ref.exons:
                    if e.start <= trec.txStart < trec.txEnd <= e.end:
                        isoform_hit.str_class = "incomplete-splice_match"
                        isoform_hit.subtype = "mono-exon"
                        isoform_hit.modify(ref.id, ref.gene, diff_tss, diff_tts, ref.length, ref.exonCount)
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
                            isoform_hit.str_class = "novel_in_catalog"
                            isoform_hit.subtype = "mono-exon_by_intron_retention"
                            
                            isoform_hit.modify("novel", ref.gene, 'NA', 'NA', ref.length, ref.exonCount)
                            get_gene_diff_tss_tts(isoform_hit,trec,start_ends_by_gene)

                            # return isoform_hit

                # if we get to here, means not ISM, so just add a ref gene and categorize further later
                isoform_hit.genes.append(ref.gene)
            if isoform_hit.str_class == "novel_in_catalog":
                return isoform_hit
    get_gene_diff_tss_tts(isoform_hit,trec,start_ends_by_gene)
    isoform_hit.genes.sort(key=lambda x: start_ends_by_gene[x]['begin'])
    return isoform_hit


def novelIsoformsKnownGenes(isoforms_hit, trec, junctions_by_chr, junctions_by_gene):
    """
    At this point: definitely not FSM or ISM, see if it is NIC, NNC, or fusion
    :return isoforms_hit: updated isoforms hit (myQueryTranscripts object)
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

    ref_genes = list(set(isoforms_hit.genes))

    # at this point, we have already found matching genes/transcripts
    # hence we do not need to update refLen or refExon
    # or tss_diff and tts_diff (always set to "NA" for non-FSM/ISM matches)
    
    isoforms_hit.transcripts = ["novel"]
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
            isoforms_hit.str_class="novel_in_catalog"
            if all_junctions_in_hit_ref:
                isoforms_hit.subtype = "combination_of_known_junctions"
            else:
                isoforms_hit.subtype = "combination_of_known_splicesites"
        else:
            isoforms_hit.str_class="novel_not_in_catalog"
            isoforms_hit.subtype = "at_least_one_novel_splicesite"
    else: # see if it is fusion
        # list of a ref junctions from all genes, including potential shared junctions
        # NOTE: some ref genes could be mono-exonic so no junctions
        all_ref_junctions = list(itertools.chain(junctions_by_gene[ref_gene] for ref_gene in ref_genes if ref_gene in junctions_by_gene))

        # (junction index) --> number of refs that have this junction
        junction_ref_hit = dict((i, all_ref_junctions.count(junc)) for i,junc in enumerate(trec.junctions))

        # if the same query junction appears in more than one of the hit references, it is not a fusion
        if max(junction_ref_hit.values()) > 1:
            isoforms_hit.str_class = "moreJunctions"
        else:
            isoforms_hit.str_class = "fusion"
            isoforms_hit.subtype = "mono-exon" if trec.exonCount==1 else "multi-exon"

    if has_intron_retention():
        isoforms_hit.subtype = "intron_retention"

    return isoforms_hit

def associationOverlapping(isoforms_hit, trec, junctions_by_chr):
    """
    Classify the given isoform based on its overlap with known genes and junctions.

    Parameters:
    isoforms_hit (object): An object representing the isoform hit, which will be updated with classification information.
    trec (object): An object representing the transcript record, containing information such as exon count, chromosome, start, and end positions.
    junctions_by_chr (dict): A dictionary where keys are chromosome names and values are dictionaries containing junction information.

    Returns:
    object: The updated isoforms_hit object with classification information added.
    """
    # at this point: definitely not FSM or ISM or NIC or NNC
    # possibly (in order of preference assignment):
    #  - antisense  (on opp strand of a known gene)
    #  - genic (overlaps a combination of exons and introns), ignore strand
    #  - genic_intron  (completely within an intron), ignore strand
    #  - intergenic (does not overlap a gene at all), ignore strand

    isoforms_hit.str_class = "intergenic"
    isoforms_hit.transcripts = ["novel"]
    isoforms_hit.subtype = "mono-exon" if trec.exonCount==1 else "multi-exon"

    #if trec.id.startswith('PB.37872.'):
    #    pdb.set_trace()
    if len(isoforms_hit.genes) == 0:
        # completely no overlap with any genes on the same strand
        # check if it is anti-sense to a known gene, otherwise it's genic_intron or intergenic
        if len(isoforms_hit.AS_genes) == 0:
            if trec.chrom in junctions_by_chr:
                # no hit even on opp strand
                # see if it is completely contained within a junction
                da_pairs = junctions_by_chr[trec.chrom]['da_tree'].find(trec.txStart, trec.txEnd)
                for junction in da_pairs:
                    if junction[0] <= trec.txStart <= trec.txEnd <= junction[1]:
                        isoforms_hit.str_class = "genic_intron"
                        break
            else:
                pass # remain intergenic
        else:
            # hits one or more genes on the opposite strand
            isoforms_hit.str_class = "antisense"
            isoforms_hit.genes = ["novelGene_{g}_AS".format(g=g) for g in isoforms_hit.AS_genes]
    else:
        # TODO: REmove this comments??
        # (Liz) used to put NNC here - now just genic
        isoforms_hit.str_class = "genic"
        # overlaps with one or more genes on the same strand
        #if trec.exonCount >= 2:
        #    # multi-exon and has a same strand gene hit, must be NNC
        #    isoforms_hit.str_class = "novel_not_in_catalog"
        #    isoforms_hit.subtype = "at_least_one_novel_splicesite"
        #else:
        #    # single exon, must be genic
        #    isoforms_hit.str_class = "genic"

    return isoforms_hit
