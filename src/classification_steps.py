
from bx.intervals import Interval
from collections import defaultdict
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

        if isoform_hit.str_class in ("anyKnownJunction", "anyKnownSpliceSite"):
            # not FSM or ISM --> see if it is NIC, NNC, or fusion
            isoform_hit = novelIsoformsKnownGenes(isoform_hit, rec, junctions_by_chr, junctions_by_gene)
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


def write_junction_info(trec, junctions_by_chr, accepted_canonical_sites, indelInfo, genome_dict,
                        fout, covInf=None, covNames=None, phyloP_reader=None):
    """
    :param trec: query isoform genePredRecord
    :param junctions_by_chr: dict of chr -> {'donors': <sorted list of donors>, 'acceptors': <sorted list of acceptors>, 'da_pairs': <sorted list of junctions>}
    :param accepted_canonical_sites: list of accepted canonical splice sites
    :param indelInfo: indels near junction information, dict of pbid --> list of junctions near indel (in Interval format)
    :param genome_dict: genome fasta dict
    :param fout: DictWriter handle
    :param covInf: (optional) junction coverage information, dict of (chrom,strand) -> (0-based start,1-based end) -> dict of {sample -> (unique, multi) read count}
    :param covNames: (optional) list of sample names for the junction coverage information
    :param phyloP_reader: (optional) dict of (chrom,0-based coord) --> phyloP score

    Write a record for each junction in query isoform
    """
    # go through each trec junction
    for junction_index, (d, a) in enumerate(trec.junctions):
        # NOTE: donor just means the start, not adjusted for strand
        # Check if the chromosome of the transcript has any annotation by the reference
        # create a list in case there are chromosomes present in the input but not in the annotation dictionary junctions_by_chr
        missing_chr=[]
        junction_cat = "novel"
        if (trec.chrom in junctions_by_chr) and (trec.chrom not in missing_chr):

            if ((d,a) in junctions_by_chr[trec.chrom]['da_pairs'][trec.strand]):
                junction_cat = "known"
                min_diff_s = min_diff_e = 0
            else:
                # Find the closest junction start site
                min_diff_s = -find_closest_in_list(junctions_by_chr[trec.chrom]['donors'], d)
                # find the closest junction end site
                min_diff_e = find_closest_in_list(junctions_by_chr[trec.chrom]['acceptors'], a)
            
        else:
            # if there is no record in the reference of junctions in this chromosome, minimum distances will be NA
            # add also new chromosome to the junctions_by_chr with one dummy SJ d=1, a=2
            if trec.chrom not in missing_chr:
                missing_chr.append(trec.chrom)
            min_diff_s = float("NaN")
            min_diff_e = float("NaN")

        splice_site = trec.get_splice_site(genome_dict, junction_index)

        indel_near_junction = "NA"
        if indelInfo is not None:
            indel_near_junction = "TRUE" if (trec.id in indelInfo and Interval(d,a) in indelInfo[trec.id]) else "FALSE"

        sample_cov = defaultdict(lambda: (0,0))  # sample -> (unique, multi) count for this junction
        if covInf is not None:
            sample_cov = covInf[(trec.chrom, trec.strand)][(d,a)]

        # if phyloP score dict exists, give the triplet score of (last base in donor exon), donor site -- similarly for acceptor
        phyloP_start, phyloP_end = 'NA', 'NA'
        if phyloP_reader is not None:
            phyloP_start = ",".join([str(x) for x in [phyloP_reader.get_pos(trec.chrom, d-1), phyloP_reader.get_pos(trec.chrom, d), phyloP_reader.get_pos(trec.chrom, d+1)]])
            phyloP_end = ",".join([str(x) for x in [phyloP_reader.get_pos(trec.chrom, a-1), phyloP_reader.get_pos(trec.chrom, a),
                                              phyloP_reader.get_pos(trec.chrom, a+1)]])

        qj = {'isoform': trec.id,
              'junction_number': "junction_"+str(junction_index+1),
              "chrom": trec.chrom,
              "strand": trec.strand,
              "genomic_start_coord": d+1,  # write out as 1-based start
              "genomic_end_coord": a,      # already is 1-based end
              "transcript_coord": "?????",  # this is where the exon ends w.r.t to id sequence, ToDo: implement later
              "junction_category": junction_cat,
              "start_site_category": "known" if min_diff_s==0 else "novel",
              "end_site_category": "known" if min_diff_e==0 else "novel",
              "diff_to_Ref_start_site": min_diff_s if min_diff_s==min_diff_s else "NA", # check if min_diff is actually nan
              "diff_to_Ref_end_site": min_diff_e if min_diff_e==min_diff_e else "NA",   # check if min_diff is actually nan
              "bite_junction": "TRUE" if ((min_diff_s<0 or min_diff_e<0) and not(min_diff_s>0 or min_diff_e>0)) else "FALSE",
              "splice_site": splice_site,
              "canonical": "canonical" if splice_site in accepted_canonical_sites else "non_canonical",
              "RTS_junction": "????", # First write ???? in _tmp, later is TRUE/FALSE
              "indel_near_junct": indel_near_junction,
              "phyloP_start": phyloP_start,
              "phyloP_end": phyloP_end,
              "sample_with_cov": sum([cov_uniq>0 for (cov_uniq,_) in sample_cov.values()]) if covInf is not None else "NA",
              "total_coverage_unique": sum([cov_uniq for (cov_uniq,_ ) in sample_cov.values()]) if covInf is not None else "NA",
              "total_coverage_multi": sum([cov_multi for (_,cov_multi ) in sample_cov.values()]) if covInf is not None else "NA"}

        if covInf is not None:
            for sample in covNames:
                cov_uniq, cov_multi = sample_cov[sample]
                qj[sample+'_unique'] = str(cov_uniq)
                qj[sample+'_multi'] = str(cov_multi)

        fout.writerow(qj)
