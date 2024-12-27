from bx.intervals import IntervalTree

def calc_overlap(s1, e1, s2, e2):
    """returns the overlap between two intervals"""
    if s1=='NA' or s2=='NA': return 0
    if s1 > s2:
        s1, e1, s2, e2 = s2, e2, s1, e1
    return max(0, min(e1,e2)-max(s1,s2) + 1) # Wouldn't this need a +1?

def calc_splicesite_agreement(query_exons, ref_exons):
    """Determines the number of splice sites that agree between query and reference"""
    q_sites = set(site for exon in query_exons for site in (exon.start, exon.end))
    return sum(1 for exon in ref_exons for site in (exon.start, exon.end) if site in q_sites)

def calc_exon_overlap(query_exons, ref_exons):
    """Determines the number of bases that overlap between query and reference"""
    query_ranges = set(b for e in query_exons for b in range(e.start, e.end))
    return sum(1 for e in ref_exons for b in range(e.start, e.end) if b in query_ranges)

def get_diff_tss_tts(trec, ref):
    """
    Calculate the differences between the Transcript Start Site (TSS) and 
    Transcript Termination Site (TTS) of two transcripts.

    For positive strand transcripts:
    - TSS is the difference between the reference's start and the transcript's start.
    - TTS is the difference between the transcript's end and the reference's end.
    Positive values indicate elongation, negative values indicate shortening.

    For negative strand transcripts:
    - The start and end are reversed (i.e., `start` corresponds to the `end` in positive strand, and vice versa).
    - TTS is calculated as the difference between the reference's start and the transcript's start.
    - TSS is the difference between the transcript's end and the reference's end.
    Again, positive values indicate elongation, negative values indicate shortening.

    Parameters:
    - trec: The transcript object representing the query transcript.
    - ref: The reference transcript object.

    Returns:
    - diff_tss (int): Difference in the transcript start sites (positive for elongation, negative for shortening).
    - diff_tts (int): Difference in the transcript termination sites (positive for elongation, negative for shortening).
    """
    if trec.strand == '+':
        diff_tss = ref.txStart - trec.txStart
        diff_tts = trec.txEnd  - ref.txEnd
    else:
        diff_tts = ref.txStart - trec.txStart
        diff_tss = trec.txEnd  - ref.txEnd
    return diff_tss, diff_tts


def get_gene_diff_tss_tts(isoform_hit,trec,start_ends_by_gene):
    """
    Calculate the nearest transcription start site (TSS) and transcription termination site (TTS) 
    differences for a given isoform hit relative to all isoforms of the gene.

    Args:
        isoform_hit: An object representing the isoform hit, which contains information about 
                        the genes it hits and their start/end sites.
        trec: The transcript record object representing the query transcript.
        start_ends_by_gene: A dictionary containing the start and end sites for all isoforms of a gene.

    Modifies:
        isoform_hit: Updates the `tss_gene_diff` and `tts_gene_diff` attributes with the nearest 
                        TSS and TTS differences, respectively. If no valid difference is found, 
                        the attribute is set to 'NA'.
    """
    # now that we know the reference (isoform) it hits
    # add the nearest start/end site for that gene (all isoforms of the gene)
    nearest_start_diff, nearest_end_diff = float('inf'), float('inf')
    for ref_gene in isoform_hit.genes:
        for x in start_ends_by_gene[ref_gene]['begin']:
            d =  x - trec.txStart
            if abs(d) < abs(nearest_start_diff):
                nearest_start_diff = d
        for x in start_ends_by_gene[ref_gene]['end']:
            d = trec.txEnd - x
            if abs(d) < abs(nearest_end_diff):
                nearest_end_diff = d

    if trec.strand == '+':
        isoform_hit.tss_gene_diff = nearest_start_diff if nearest_start_diff!=float('inf') else 'NA'
        isoform_hit.tts_gene_diff = nearest_end_diff if nearest_end_diff!=float('inf') else 'NA'
    else:
        isoform_hit.tss_gene_diff = nearest_end_diff if nearest_start_diff!=float('inf') else 'NA'
        isoform_hit.tts_gene_diff = nearest_start_diff if nearest_end_diff!=float('inf') else 'NA'

def categorize_incomplete_matches(trec, ref):
    """
    intron_retention --- at least one trec exon covers at least two adjacent ref exons
    complete --- all junctions agree and is not IR
    5prime_fragment --- all junctions agree but trec has less 5' exons. The isoform is a 5' fragment of the reference transcript
    3prime_fragment --- all junctions agree but trec has less 3' exons. The isoform is a 3' fragment of the reference transcript
    internal_fragment --- all junctions agree but trec has less 5' and 3' exons
    """
    # check intron retention
    ref_exon_tree = IntervalTree()
    for i,e in enumerate(ref.exons): ref_exon_tree.insert(e.start, e.end, i)
    for e in trec.exons:
        if len(ref_exon_tree.find(e.start, e.end)) > 1: # multiple ref exons covered
            return "intron_retention"

    agree_front = trec.junctions[0]==ref.junctions[0]
    agree_end   = trec.junctions[-1]==ref.junctions[-1]
    if agree_front:
        if agree_end:
            return "complete"
        else: # front agrees, end does not
            return ("5prime_fragment" if trec.strand=='+' else '3prime_fragment')
    else:
        if agree_end: # front does not agree, end agrees
            return ("3prime_fragment" if trec.strand=='+' else '5prime_fragment')
        else: # TODO: check if the case of the test would be possible
            return "internal_fragment"
        
def full_splice_match_subtype(diff_tss,diff_tts):
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
    return subtype