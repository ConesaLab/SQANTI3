
from bx.intervals import Interval
from collections import defaultdict
from src.classification_utils import add_coding_info
from src.qc_classes import myQueryProteins
from src.utils import find_closest_in_list, find_polyA_motif
from src.config import seqid_fusion
from src.classification_classifiers import (
    transcriptsKnownSpliceSites, novelIsoformsKnownGenes, associationOverlapping
)

def classify_isoform(rec, refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr, junctions_by_gene,
                     start_ends_by_gene, genome_dict, isoform_hits_name=None,window=20):
        # Find best reference hit
        isoform_hit = transcriptsKnownSpliceSites(isoform_hits_name, refs_1exon_by_chr, refs_exons_by_chr, 
                                                  start_ends_by_gene, rec, genome_dict, nPolyA=window)

        if isoform_hit.structural_category in ("anyKnownJunction", "anyKnownSpliceSite"):
            # not FSM or ISM --> see if it is NIC, NNC, or fusion
            isoform_hit = novelIsoformsKnownGenes(isoform_hit, rec, junctions_by_chr, junctions_by_gene)
        elif isoform_hit.structural_category in ("", "antisense","geneOverlap"):
            # possibly NNC, genic, genic intron, anti-sense, or intergenic
            isoform_hit = associationOverlapping(isoform_hit, rec, junctions_by_chr)

        return isoform_hit

def process_cage_peak_info(isoform_hit,rec, cage_peak_obj):
    if rec.strand == '+':
        within_CAGE, dist_CAGE = cage_peak_obj.find(rec.chrom, rec.strand, rec.txStart)
    else:
        within_CAGE, dist_CAGE = cage_peak_obj.find(rec.chrom, rec.strand, rec.txEnd)
    isoform_hit.within_CAGE_peak = within_CAGE
    isoform_hit.dist_to_CAGE_peak = dist_CAGE


def process_polya_peak_info(isoform_hit,rec, polya_peak_obj):
    if rec.strand == '+':
        within_polyA_site, dist_polyA_site = polya_peak_obj.find(rec.chrom, rec.strand, rec.txEnd)
    else:
        within_polyA_site, dist_polyA_site = polya_peak_obj.find(rec.chrom, rec.strand, rec.txStart)
    isoform_hit.within_polyA_site = within_polyA_site
    isoform_hit.dist_to_polyA_site = dist_polyA_site

def find_polya_motif_info(isoform_hit, rec, genome_dict, polyA_motif_list):
    if rec.strand == '+':
        polyA_motif, polyA_dist, polyA_motif_found = find_polyA_motif(str(genome_dict[rec.chrom][rec.txEnd-50:rec.txEnd].seq), polyA_motif_list)
    else:
        polyA_motif, polyA_dist, polyA_motif_found = find_polyA_motif(str(genome_dict[rec.chrom][rec.txStart:rec.txStart+50].reverse_complement().seq), polyA_motif_list)
    isoform_hit.polyA_motif = polyA_motif
    isoform_hit.polyA_dist = polyA_dist
    isoform_hit.polyA_motif_found = polyA_motif_found

def fill_cds_info(isoform_hit, rec, cdsDict, is_fusion, fusion_components):
    if is_fusion:
        fusion_gene = 'PBfusion.' + str(seqid_fusion.match(rec.id).group(1))
        if fusion_gene not in cdsDict:
            return
        
        cds_info = cdsDict[fusion_gene]
        rec_start, rec_end = fusion_components[rec.id]
        rec_len = rec_end - rec_start + 1
        orf_start, orf_end = cds_info.cds_start, cds_info.cds_end

        # CASE 1: component starts inside ORF
        if orf_start <= rec_start < orf_end:
            cds_start = 1
            cds_end = min(rec_len, orf_end - rec_start + 1)
            protein_len = (cds_end - cds_start) // 3
            offset = (rec_start - orf_start) // 3
            seq_end = min(offset + protein_len, len(cds_info.protein_seq))
            protein_seq = cds_info.protein_seq[offset:seq_end]

        # CASE 2: ORF starts inside component
        elif rec_start <= orf_start < rec_end:
            cds_start = orf_start - rec_start
            cds_end = (rec_end - rec_start + 1
                       if orf_end >= rec_end else orf_end - rec_start + 1)
            protein_len = (cds_end - cds_start) // 3
            seq_end = min(protein_len, len(cds_info.protein_seq))
            protein_seq = cds_info.protein_seq[:seq_end]

        else:
            return  # No overlap — leave isoform_hit unchanged

        # Create a temporary myQueryProteins object with the mapped info
        fusion_cds = myQueryProteins(
            cds_start=cds_start,
            cds_end=cds_end,
            protein_length=protein_len,
            protein_seq=protein_seq,
            proteinID=fusion_gene,
            psauron_score=cds_info.psauron_score,
            cds_type=cds_info.cds_type
        )
        add_coding_info(isoform_hit, fusion_cds)

    elif rec.id in cdsDict:
        add_coding_info(isoform_hit, cdsDict[rec.id])


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
    #cdsDict[rec.id].cds_genomic_start = m[cdsDict[rec.id].cds_start-1] + 1  # make it 1-based
    #cdsDict[rec.id].cds_genomic_end   = m[cdsDict[rec.id].cds_end-1] + 1    # make it 1-based

def detect_nmd(isoform_hit, rec):
    # NMD detection
    # if + strand, see if CDS stop is before the last junction
    if len(rec.junctions) > 0:
        if rec.strand == '+':
            dist_to_last_junc = isoform_hit.CDS_genomic_end - rec.junctions[-1][0]
        else: # - strand
            dist_to_last_junc = rec.junctions[0][1] - isoform_hit.CDS_genomic_end
        isoform_hit.predicted_NMD = "TRUE" if dist_to_last_junc < -50 else "FALSE"

