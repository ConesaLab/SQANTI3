import sys
from csv import DictWriter

from src.classification_steps import (
    assign_genomic_coordinates, classify_isoform, 
    detect_nmd, fill_orf_info, find_polya_motif_info,
    process_cage_peak_info, process_polya_peak_info,
    write_junction_info
) # type: ignore
from src.utils import alphanum_key
from src.config import  FIELDS_CLASS
from src.module_logging import qc_logger


def isoform_classification_pipeline(
        sites,window,is_fusion, isoforms_by_chr, refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr,
        junctions_by_gene, start_ends_by_gene, genome_dict, indelsJunc, orfDict,
        outputClassPath, outputJuncPath, fusion_components,isoform_hits_name,SJcovNames, 
        SJcovInfo, fields_junc_cur,ratio_TSS_dict, cage_peak_obj, polya_peak_obj,
        polyA_motif_list, phyloP_reader
        ):
    # running classification
    qc_logger.info("**** Performing Classification of Isoforms")

    accepted_canonical_sites = list(sites.split(","))
    # Creates a temporary file to write the classification and junction results
    handle_class = open(outputClassPath+"_tmp", "w")
    fout_class = DictWriter(handle_class, fieldnames=FIELDS_CLASS, delimiter='\t')
    fout_class.writeheader()

    handle_junc = open(outputJuncPath+"_tmp", "w")
    fout_junc = DictWriter(handle_junc, fieldnames=fields_junc_cur, delimiter='\t')
    fout_junc.writeheader()

    isoforms_info = {}
    for _,records in isoforms_by_chr.items():
        for rec in records:
            # Find best reference hit
            isoform_hit = classify_isoform(rec, refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr,
                                            junctions_by_gene, start_ends_by_gene, genome_dict, isoform_hits_name,window)
            # write out junction information
            write_junction_info(rec, junctions_by_chr, accepted_canonical_sites, indelsJunc, 
                                genome_dict, fout_junc, SJcovInfo, SJcovNames, phyloP_reader)

            # look at Cage Peak info (if available)
            if cage_peak_obj is not None:
                process_cage_peak_info(isoform_hit, rec, cage_peak_obj)
            
            # look at PolyA Peak info (if available)
            if polya_peak_obj is not None:
                process_polya_peak_info(isoform_hit, rec, polya_peak_obj)
            
            # polyA motif finding: look within 50 bp upstream of 3' end for the highest ranking polyA motif signal (user provided)
            if polyA_motif_list is not None:
                find_polya_motif_info(isoform_hit, rec, genome_dict, polyA_motif_list)

            # Fill in ORF/coding info and NMD detection
            fill_orf_info(isoform_hit, rec, orfDict, is_fusion, fusion_components)

            # Assign the genomic coordinates of the CDS start and end
            if isoform_hit.coding == "coding":
                assign_genomic_coordinates(isoform_hit, rec)

            if isoform_hit.CDS_genomic_end!='NA':
                detect_nmd(isoform_hit, rec)
            isoforms_info[rec.id] = isoform_hit
            fout_class.writerow(isoform_hit.as_dict())

    # Sort isoforms_info by the chromosome and then id.
    # Take into account that the ID is a string with numbers

    iso_keys = sorted(isoforms_info.keys(), key=lambda x: (alphanum_key(isoforms_info[x].chrom),
                                                           alphanum_key(isoforms_info[x].id)))
    isoforms_info = {key: isoforms_info[key] for key in iso_keys}
    
    handle_class.close()
    handle_junc.close()
    return (isoforms_info, ratio_TSS_dict)
