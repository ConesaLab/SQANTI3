
from .utils import find_polyA_motif
from .config import seqid_fusion
from .classification_classifiers import (
    transcriptsKnownSpliceSites, novelIsoformsKnownGenes, associationOverlapping
)

def classify_isoform(rec, refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr, junctions_by_gene,
                     start_ends_by_gene, genome_dict, isoform_hits_name=None,window=20):
        # Find best reference hit
        isoform_hit = transcriptsKnownSpliceSites(isoform_hits_name, refs_1exon_by_chr, refs_exons_by_chr, 
                                                  start_ends_by_gene, rec, genome_dict, nPolyA=window)

        if isoform_hit.str_class in ("anyKnownJunction", "anyKnownSpliceSite"):
            # not FSM or ISM --> see if it is NIC, NNC, or fusion
            isoform_hit = novelIsoformsKnownGenes(isoform_hit, rec, junctions_by_chr, junctions_by_gene, 
                                                  start_ends_by_gene)
        elif isoform_hit.str_class in ("", "geneOverlap"):
            # possibly NNC, genic, genic intron, anti-sense, or intergenic
            isoform_hit = associationOverlapping(isoform_hit, rec, junctions_by_chr)

        return isoform_hit

def process_cage_peak_info(isoform_hit,rec, cage_peak_obj):
    if rec.strand == '+':
        within_CAGE, dist_CAGE = cage_peak_obj.find(rec.chrom, rec.strand, rec.txStart)
    else:
        within_CAGE, dist_CAGE = cage_peak_obj.find(rec.chrom, rec.strand, rec.txEnd)
    isoform_hit.within_CAGE = within_CAGE
    isoform_hit.dist_CAGE = dist_CAGE


def process_polya_peak_info(isoform_hit,rec, polya_peak_obj):
    if rec.strand == '+':
        within_polyA_site, dist_polyA_site = polya_peak_obj.find(rec.chrom, rec.strand, rec.txEnd)
    else:
        within_polyA_site, dist_polyA_site = polya_peak_obj.find(rec.chrom, rec.strand, rec.txStart)
    isoform_hit.within_polyA_site = within_polyA_site
    isoform_hit.dist_polyA_site = dist_polyA_site

def find_polya_motif_info(isoform_hit, rec, genome_dict, polyA_motif_list):
    if rec.strand == '+':
        polyA_motif, polyA_dist, polyA_motif_found = find_polyA_motif(str(genome_dict[rec.chrom][rec.txEnd-50:rec.txEnd].seq), polyA_motif_list)
    else:
        polyA_motif, polyA_dist, polyA_motif_found = find_polyA_motif(str(genome_dict[rec.chrom][rec.txStart:rec.txStart+50].reverse_complement().seq), polyA_motif_list)
    isoform_hit.polyA_motif = polyA_motif
    isoform_hit.polyA_dist = polyA_dist
    isoform_hit.polyA_motif_found = polyA_motif_found

def fill_orf_info(isoform_hit, rec, orfDict, is_fusion, fusion_components):
    if is_fusion:
        # fusion - special case handling, need to see which part of the ORF this segment falls on
        fusion_gene = 'PBfusion.' + str(seqid_fusion.match(rec.id).group(1))
        rec_component_start, rec_component_end = fusion_components[rec.id]
        rec_len = rec_component_end - rec_component_start + 1
        if fusion_gene in orfDict:
            orf_start, orf_end = orfDict[fusion_gene].cds_start, orfDict[fusion_gene].cds_end
            if orf_start <= rec_component_start < orf_end:
                isoform_hit.CDS_start = 1
                isoform_hit.CDS_end = min(rec_len, orf_end - rec_component_start + 1)
                isoform_hit.ORFlen = (isoform_hit.CDS_end - isoform_hit.CDS_start) / 3
                _s = (rec_component_start - orf_start) // 3
                _e = min(int(_s + isoform_hit.ORFlen), len(orfDict[fusion_gene].orf_seq))
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
    elif rec.id in orfDict: # this will never be true for fusion, so the above code seg runs instead
        isoform_hit.coding = "coding"
        isoform_hit.ORFlen = orfDict[rec.id].orf_length
        isoform_hit.CDS_start = orfDict[rec.id].cds_start
        isoform_hit.CDS_end = orfDict[rec.id].cds_end
        isoform_hit.ORFseq = orfDict[rec.id].orf_seq

def assign_genomic_coordinates(isoform_hit, rec):
    m = {}
    if rec.strand == '+':
        i = 0
        for exon in rec.exons:
            for c in range(exon.start, exon.end):
                m[i] = c
                i += 1
    else:
        i = 0
        for exon in rec.exons:
            for c in range(exon.start, exon.end):
                m[rec.length - i - 1] = c
                i += 1

    isoform_hit.CDS_genomic_start = m[isoform_hit.CDS_start-1] + 1  # make it 1-based
    # NOTE: if using --orf_input, it is possible to see discrepancy between the exon structure
    # provided by GFF and the input ORF. For now, just shorten it
    isoform_hit.CDS_genomic_end = m[min(isoform_hit.CDS_end-1, max(m))] + 1    # make it 1-based
    #orfDict[rec.id].cds_genomic_start = m[orfDict[rec.id].cds_start-1] + 1  # make it 1-based
    #orfDict[rec.id].cds_genomic_end   = m[orfDict[rec.id].cds_end-1] + 1    # make it 1-based

def detect_nmd(isoform_hit, rec):
    # NMD detection
    # if + strand, see if CDS stop is before the last junction
    if len(rec.junctions) > 0:
        if rec.strand == '+':
            dist_to_last_junc = isoform_hit.CDS_genomic_end - rec.junctions[-1][0]
        else: # - strand
            dist_to_last_junc = rec.junctions[0][1] - isoform_hit.CDS_genomic_end
        isoform_hit.is_NMD = "TRUE" if dist_to_last_junc < -50 else "FALSE"
