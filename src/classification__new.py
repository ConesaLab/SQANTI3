import os
import sys
import csv
import bisect
import itertools
from csv import DictWriter
from collections import defaultdict, namedtuple
from bx.intervals import IntervalTree # type: ignore

from .utilities.cupcake.tofu.compare_junctions import compare_junctions
from .utilities.cupcake.sequence.BED import LazyBEDPointReader
from .utilities.short_reads import get_bam_header,get_ratio_TSS, get_TSS_bed

from .qc_classes import myQueryTranscripts, CAGEPeak, PolyAPeak
from .helpers import write_junctionInfo, get_isoform_hits_name
from .config import FIELDS_JUNC, FIELDS_CLASS, seqid_fusion
from .parsers import get_fusion_component
from .utils import find_polyA_motif
from .classification_utils import (
    SJ_coverage, TSS_ratio_calculation, calc_exon_overlap, calc_splicesite_agreement, full_splice_match_subtype, get_diff_tss_tts,
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
            main_type, subtype = "", "no_subcategory"
            diff_tss, diff_tts = "NA", "NA"
            splicesite_agreement = None
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
                exon_overlap = calc_exon_overlap(trec.exons, ref.exons)
                
                if ref.exonCount == 1: # mono-exonic reference, handle specially here
                    # TODO: calc exon overlap is called more than twice. Create a variable to store the value 
                    if exon_overlap > 0 and cat_ranking[main_type] < cat_ranking["geneOverlap"]: #CHECK wouldn't it be better with just a number?
                        main_type = "geneOverlap"
                        subtype = "mono-exon"
                        splicesite_agreement = 0
                        diff_tss, diff_tts = "NA", "NA"


                else: # multi-exonic reference
                    match_type = compare_junctions(trec, ref, internal_fuzzy_max_dist=0, max_5_diff=999999, max_3_diff=999999)
                    #TODO: Error handling in logs
                    if match_type not in ('exact', 'subset', 'partial', 'concordant', 'super', 'nomatch'):
                        raise Exception("Unknown match category {0}!".format(match_type))

                    diff_tss, diff_tts = get_diff_tss_tts(trec, ref)
                    splicesite_agreement = calc_splicesite_agreement(trec.exons, ref.exons)

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
                        if cat_ranking[main_type] < cat_ranking["full-splice_match"] or \
                                                    abs(diff_tss)+abs(diff_tts) < isoform_hit.get_total_diff(): # FIX: check this number
                            main_type = "full-splice_match"
                            subtype = full_splice_match_subtype(diff_tss,diff_tts)

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
                        if cat_ranking[main_type] < cat_ranking["incomplete-splice_match"] or \
                            (main_type=='incomplete-splice_match' and abs(diff_tss)+abs(diff_tts) < isoform_hit.get_total_diff()):
                            main_type = "incomplete-splice_match"
   
                    # #######################################################
                    # Some kind of junction match that isn't ISM/FSM
                    # #######################################################
                    elif match_type in ('partial', 'concordant', 'super'):
                        q_exon_d = abs(trec.exonCount - ref.exonCount)
                        main_type = "anyKnownJunction"
                        if cat_ranking[main_type] < cat_ranking["anyKnownJunction"] or \
                                (main_type=='anyKnownJunction' and splicesite_agreement  > isoform_hit.q_splicesite_hit) or \
                                (main_type=='anyKnownJunction' and splicesite_agreement == isoform_hit.q_splicesite_hit and exon_overlap > isoform_hit.q_exon_overlap) or \
                                (isoform_hit.str_class=='anyKnownJunction' and splicesite_agreement == isoform_hit.q_splicesite_hit and q_exon_d < abs(trec.exonCount-isoform_hit.refExons)):
                            subtype = "no_subcategory"

                    else: # must be nomatch
                        assert match_type == 'nomatch'
                        # at this point, no junction overlap, but may be a single splice site (donor or acceptor) match?
                        # also possibly just exonic (no splice site) overlap
                        if cat_ranking[main_type] < cat_ranking["anyKnownSpliceSite"] and calc_splicesite_agreement(trec.exons, ref.exons) > 0:
                            main_type = "anyKnownSpliceSite"
                            subtype = "no_subcategory"
                            diff_tss, diff_tts = "NA", "NA"

                        if main_type == "": # still not hit yet, check exonic overlap
                            if cat_ranking[main_type] < cat_ranking["geneOverlap"] and calc_exon_overlap(trec.exons, ref.exons) > 0:
                                main_type = "geneOverlap"
                                subtype = "no_subcategory"
                                diff_tss, diff_tts = "NA", "NA"

                try:
                    isoform_hit = myQueryTranscripts(trec.id, diff_tss, diff_tts, trec.exonCount, trec.length,
                                                str_class=main_type,
                                                subtype=subtype,
                                                chrom=trec.chrom,
                                                strand=trec.strand,
                                                genes=[ref.gene],
                                                transcripts=[ref.id],
                                                refLen=ref.length,
                                                refExons=ref.exonCount,
                                                refStart=ref.txStart,
                                                refEnd=ref.txEnd,
                                                q_splicesite_hit=splicesite_agreement,
                                                q_exon_overlap=exon_overlap,
                                                percAdownTTS=str(percA),
                                                seqAdownTTS=seq_downTTS)
                except Exception as e:
                    print(trec.id,match_type)
                    print(f"Error: {e}")
                    sys.exits(1)


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


def novelIsoformsKnownGenes(isoforms_hit, trec, junctions_by_chr, junctions_by_gene, start_ends_by_gene):
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
        for e in trec.exons:
            m = bisect.bisect_left(junctions_by_chr[trec.chrom]['da_pairs'], (e.start, e.end))
            if m < len(junctions_by_chr[trec.chrom]['da_pairs']) and e.start <= junctions_by_chr[trec.chrom]['da_pairs'][m][0] < junctions_by_chr[trec.chrom]['da_pairs'][m][1] < e.end:
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
                da_pairs = junctions_by_chr[trec.chrom]['da_pairs']
                i = bisect.bisect_left(da_pairs, (trec.txStart, trec.txEnd))
                for j in range(i-1, min(i+1, len(da_pairs)-1)):
                    if da_pairs[j][0] <= trec.txStart <= trec.txEnd <= da_pairs[j][1]:
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


def preprocess_isoform_data(args, corrGTF):
    """
    Preprocesses isoform data and initializes various objects for isoform classification.

    Parameters:
    args (argparse.Namespace): Command-line arguments.
    corrGTF (str): Path to the corrected GTF file.

    Returns:
    tuple: A tuple containing:
        - isoform_hits_name (str): Name of the isoform hits file.
        - star_out (str): Path to the STAR output directory.
        - star_index (str): Path to the STAR index directory.
        - SJcovNames (list): List of splice junction coverage file names.
        - SJcovInfo (dict): Dictionary containing splice junction coverage information.
        - fields_junc_cur (list): List of fields for junction information.
        - ratio_TSS_dict (dict): Dictionary containing TSS ratio information.
        - cage_peak_obj (CAGEPeak): CAGE peak object.
        - polya_peak_obj (PolyAPeak): PolyA peak object.
        - polyA_motif_list (list): List of PolyA motifs.
        - phyloP_reader (LazyBEDPointReader): PhyloP BED reader object.
    """
    global isoform_hits_name
    isoform_hits_name = None
    fusion_components = {}
    cage_peak_obj = None
    polya_peak_obj = None
    polyA_motif_list = None
    phyloP_reader = None

    if args.is_fusion: # read GTF to get fusion components
        # ex: PBfusion.1.1 --> (1-based start, 1-based end) of where the fusion component is w.r.t to entire fusion
        fusion_components = get_fusion_component(args.isoforms)

    # TODO: Create a new module to do all of the pre-classification steps.
    # If the isoform hits are present:
    if args.isoform_hits:
        isoform_hits_name = get_isoform_hits_name(args.dir,args.output)
        with open(isoform_hits_name+'_tmp', 'w') as out_file:
            tsv_writer = csv.writer(out_file, delimiter='\t')
            tsv_writer.writerow(['Isoform', 'Isoform_length', 'Isoform_exon_number', 'Hit', 'Hit_length',
                                 'Hit_exon_number', 'Match', 'Diff_to_TSS', 'Diff_to_TTS', 'Matching_type'])

    ## read coverage files if provided
    star_out, star_index, \
        SJcovNames, SJcovInfo, fields_junc_cur \
            = SJ_coverage(args.short_reads, args.coverage, 
                          args.genome, args.dir, args.cpus)

    ## TSS ratio calculation
    ratio_TSS_dict = TSS_ratio_calculation(args.SR_bam,args.short_reads,
                                           star_out,star_index,corrGTF,args.ratio_TSS_metric)
    # CAGE peaks
    if args.CAGE_peak is not None:
        print("**** Reading CAGE Peak data.", file=sys.stdout)
        cage_peak_obj = CAGEPeak(args.CAGE_peak)

    # PolyA peaks
    if args.polyA_peak is not None:
        print("**** Reading polyA Peak data.", file=sys.stdout)
        polya_peak_obj = PolyAPeak(args.polyA_peak)

    # PolyA motif list
    if args.polyA_motif_list is not None:
        print("**** Reading PolyA motif list.", file=sys.stdout)
        polyA_motif_list = []
        for line in open(args.polyA_motif_list):
            x = line.strip().upper().replace('U', 'A')
            if any(s not in ('A','T','C','G') for s in x):
                print("PolyA motif must be A/T/C/G only! Saw: {0}. Abort!".format(x), file=sys.stderr)
                sys.exit(1)
            polyA_motif_list.append(x)

    # PhyloP
    if args.phyloP_bed is not None:
        print("**** Reading PhyloP BED file.", file=sys.stdout)
        phyloP_reader = LazyBEDPointReader(args.phyloP_bed)

    return (fusion_components,isoform_hits_name, SJcovNames, SJcovInfo, fields_junc_cur,
            ratio_TSS_dict, cage_peak_obj, polya_peak_obj, polyA_motif_list, phyloP_reader)


def isoformClassification(args, isoforms_by_chr, refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr, 
                          junctions_by_gene, start_ends_by_gene, genome_dict, indelsJunc, orfDict,
                          outputClassPath, outputJuncPath, fusion_components,isoform_hits_name,SJcovNames, 
                          SJcovInfo, fields_junc_cur,ratio_TSS_dict, cage_peak_obj, polya_peak_obj,
                          polyA_motif_list, phyloP_reader):
    # running classification
    print("**** Performing Classification of Isoforms....", file=sys.stdout)


    accepted_canonical_sites = list(args.sites.split(","))
    # Creates a temporary file to write the classification and junction results
    handle_class = open(outputClassPath+"_tmp", "w")
    fout_class = DictWriter(handle_class, fieldnames=FIELDS_CLASS, delimiter='\t')
    fout_class.writeheader()

    handle_junc = open(outputJuncPath+"_tmp", "w")
    fout_junc = DictWriter(handle_junc, fieldnames=fields_junc_cur, delimiter='\t')
    fout_junc.writeheader()

    isoforms_info = {}
    novel_gene_index = 1

    for _,records in isoforms_by_chr.items():
        for rec in records:
            # Find best reference hit
            isoform_hit = transcriptsKnownSpliceSites(isoform_hits_name, refs_1exon_by_chr, refs_exons_by_chr, start_ends_by_gene, rec, genome_dict, nPolyA=args.window)

            if isoform_hit.str_class in ("anyKnownJunction", "anyKnownSpliceSite"):
                # not FSM or ISM --> see if it is NIC, NNC, or fusion
                isoform_hit = novelIsoformsKnownGenes(isoform_hit, rec, junctions_by_chr, junctions_by_gene, start_ends_by_gene)
            elif isoform_hit.str_class in ("", "geneOverlap"):
                # possibly NNC, genic, genic intron, anti-sense, or intergenic
                isoform_hit = associationOverlapping(isoform_hit, rec, junctions_by_chr)

            # write out junction information
            write_junctionInfo(rec, junctions_by_chr, accepted_canonical_sites, indelsJunc, genome_dict, fout_junc, covInf=SJcovInfo, covNames=SJcovNames, phyloP_reader=phyloP_reader)

            if isoform_hit.str_class in ("intergenic", "genic_intron"):
                # Liz: I don't find it necessary to cluster these novel genes. They should already be always non-overlapping.
                if args.novel_gene_prefix is not None:  # used by splits to not have redundant novelGene IDs
                    isoform_hit.genes = ['novelGene_' + str(args.novel_gene_prefix) + '_' + str(novel_gene_index)]
                else:
                    isoform_hit.genes = ['novelGene_' + str(novel_gene_index)]
                isoform_hit.transcripts = ['novel']
                novel_gene_index += 1

            # look at Cage Peak info (if available)
            if cage_peak_obj is not None:
                if rec.strand == '+':
                    within_CAGE, dist_CAGE = cage_peak_obj.find(rec.chrom, rec.strand, rec.txStart)
                else:
                    within_CAGE, dist_CAGE = cage_peak_obj.find(rec.chrom, rec.strand, rec.txEnd)
                isoform_hit.within_CAGE = within_CAGE
                isoform_hit.dist_CAGE = dist_CAGE

            # look at PolyA Peak info (if available)
            if polya_peak_obj is not None:
                if rec.strand == '+':
                    within_polyA_site, dist_polyA_site = polya_peak_obj.find(rec.chrom, rec.strand, rec.txEnd)
                else:
                    within_polyA_site, dist_polyA_site = polya_peak_obj.find(rec.chrom, rec.strand, rec.txStart)
                isoform_hit.within_polyA_site = within_polyA_site
                isoform_hit.dist_polyA_site = dist_polyA_site

            # polyA motif finding: look within 50 bp upstream of 3' end for the highest ranking polyA motif signal (user provided)
            if polyA_motif_list is not None:
                if rec.strand == '+':
                    polyA_motif, polyA_dist, polyA_motif_found = find_polyA_motif(str(genome_dict[rec.chrom][rec.txEnd-50:rec.txEnd].seq), polyA_motif_list)
                else:
                    polyA_motif, polyA_dist, polyA_motif_found = find_polyA_motif(str(genome_dict[rec.chrom][rec.txStart:rec.txStart+50].reverse_complement().seq), polyA_motif_list)
                isoform_hit.polyA_motif = polyA_motif
                isoform_hit.polyA_dist = polyA_dist
                isoform_hit.polyA_motif_found = polyA_motif_found

            # Fill in ORF/coding info and NMD detection
            if args.is_fusion:
                #pdb.set_trace()
                # fusion - special case handling, need to see which part of the ORF this segment falls on
                fusion_gene = 'PBfusion.' + str(seqid_fusion.match(rec.id).group(1))
                rec_component_start, rec_component_end = fusion_components[rec.id]
                rec_len = rec_component_end - rec_component_start + 1
                if fusion_gene in orfDict:
                    orf_start, orf_end = orfDict[fusion_gene].cds_start, orfDict[fusion_gene].cds_end
                    if orf_start <= rec_component_start < orf_end:
                        isoform_hit.CDS_start = 1
                        isoform_hit.CDS_end = min(rec_len, orf_end - rec_component_start + 1)
                        isoform_hit.ORFlen = (isoform_hit.CDS_end - isoform_hit.CDS_start)/3
                        _s = (rec_component_start-orf_start)//3
                        _e = min(int(_s+isoform_hit.ORFlen), len(orfDict[fusion_gene].orf_seq))
                        isoform_hit.ORFseq = orfDict[fusion_gene].orf_seq[_s:_e]
                        isoform_hit.coding = "coding"
                    elif rec_component_start <= orf_start < rec_component_end:
                        isoform_hit.CDS_start = orf_start - rec_component_start
                        if orf_end >= rec_component_end:
                            isoform_hit.CDS_end = rec_component_end - rec_component_start + 1
                        else:
                            isoform_hit.CDS_end = orf_end - rec_component_start + 1
                        isoform_hit.ORFlen = (isoform_hit.CDS_end - isoform_hit.CDS_start) / 3
                        _e = min(int(isoform_hit.ORFlen), len(orfDict[fusion_gene].orf_seq))
                        isoform_hit.ORFseq = orfDict[fusion_gene].orf_seq[:_e]
                        isoform_hit.coding = "coding"
            elif rec.id in orfDict:  # this will never be true for fusion, so the above code seg runs instead
                isoform_hit.coding = "coding"
                isoform_hit.ORFlen = orfDict[rec.id].orf_length
                isoform_hit.CDS_start = orfDict[rec.id].cds_start  # 1-based start
                isoform_hit.CDS_end = orfDict[rec.id].cds_end      # 1-based end
                isoform_hit.ORFseq  = orfDict[rec.id].orf_seq

            # Assign the genomic coordinates of the CDS start and end
            if isoform_hit.coding == "coding":
                m = {} # transcript coord (0-based) --> genomic coord (0-based)
                if rec.strand == '+':
                    i = 0
                    for exon in rec.exons:
                        for c in range(exon.start, exon.end):
                            m[i] = c
                            i += 1
                else: # - strand
                    i = 0
                    for exon in rec.exons:
                        for c in range(exon.start, exon.end):
                            m[rec.length-i-1] = c
                            i += 1

                isoform_hit.CDS_genomic_start = m[isoform_hit.CDS_start-1] + 1  # make it 1-based
                # NOTE: if using --orf_input, it is possible to see discrepancy between the exon structure
                # provided by GFF and the input ORF. For now, just shorten it
                isoform_hit.CDS_genomic_end = m[min(isoform_hit.CDS_end-1, max(m))] + 1    # make it 1-based
                #orfDict[rec.id].cds_genomic_start = m[orfDict[rec.id].cds_start-1] + 1  # make it 1-based
                #orfDict[rec.id].cds_genomic_end   = m[orfDict[rec.id].cds_end-1] + 1    # make it 1-based


            if isoform_hit.CDS_genomic_end!='NA':
                # NMD detection
                # if + strand, see if CDS stop is before the last junction
                if len(rec.junctions) > 0:
                    if rec.strand == '+':
                        dist_to_last_junc = isoform_hit.CDS_genomic_end - rec.junctions[-1][0]
                    else: # - strand
                        dist_to_last_junc = rec.junctions[0][1] - isoform_hit.CDS_genomic_end
                    isoform_hit.is_NMD = "TRUE" if dist_to_last_junc < -50 else "FALSE"

            isoforms_info[rec.id] = isoform_hit
            fout_class.writerow(isoform_hit.as_dict())

    handle_class.close()
    handle_junc.close()
    return (isoforms_info, ratio_TSS_dict)
