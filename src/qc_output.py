from collections import defaultdict
from csv import DictReader, DictWriter, writer
import os, sys
import pickle

from bx.intervals import Interval

from src.utilities.cupcake.io.GFF import collapseGFFReader, write_collapseGFF_format
from src.commands import (
    RSCRIPTPATH, RSCRIPT_QC_REPORT, RSCRIPT_TUSCO_REPORT, run_command

)
from src.config import utilitiesPath
from src.helpers import get_isoform_hits_name, get_omitted_name
from src.module_logging import qc_logger
from src.utils import find_closest_in_list

def write_omitted_isoforms(isoforms_info, outdir,prefix,min_ref_len,is_fusion, fields_class_cur):
    if min_ref_len > 0 and not is_fusion:
        omitted_name = get_omitted_name(outdir,prefix)
        omitted_iso = {key: isoforms_info[key] for key in isoforms_info 
                       if not isoforms_info[key].refLen == 'NA' and int(isoforms_info[key].refLen) <= int(min_ref_len)}
        
        for key in omitted_iso:
            del isoforms_info[key]
        
        omitted_keys = sorted(omitted_iso.keys(), key=lambda x: (omitted_iso[x].chrom, omitted_iso[x].id))
        
        with open(omitted_name, 'w') as h:
            fout_omitted = DictWriter(h, fieldnames=fields_class_cur, delimiter='\t')
            fout_omitted.writeheader()
            for key in omitted_keys:
                fout_omitted.writerow(omitted_iso[key].as_dict())
    
    return isoforms_info

def write_classification_output(isoforms_info, outputClassPath, fields_class_cur):
    with open(outputClassPath, 'w') as h:
        header = next(iter(isoforms_info.values())).as_dict().keys() 
        fout_class = DictWriter(h, fieldnames=list(header), delimiter='\t')
        fout_class.writeheader()
        for iso_key in isoforms_info.keys():
            fout_class.writerow(isoforms_info[iso_key].as_dict())

def write_junction_output(outputJuncPath, RTS_info, fields_junc_cur):
    with open(outputJuncPath, 'w') as h:
        fout_junc = DictWriter(h, fieldnames=fields_junc_cur, delimiter='\t')
        fout_junc.writeheader()
        for r in DictReader(open(outputJuncPath+"_tmp"), delimiter='\t'):
            if r['isoform'] in RTS_info:
                r['RTS_junction'] = 'TRUE' if r['junction_number'] in RTS_info[r['isoform']] else 'FALSE'
            fout_junc.writerow(r)

def write_isoform_hits(outdir,prefix, isoforms_info):
    fields_hits = ['Isoform', 'Isoform_length', 'Isoform_exon_number', 'Hit', 'Hit_length', 'Hit_exon_number', 'Match', 'Diff_to_TSS', 'Diff_to_TTS', 'Matching_type']
    isoform_hits_name = get_isoform_hits_name(outdir, prefix)
    with open(isoform_hits_name, 'w') as h:
        fout_hits = DictWriter(h, fieldnames=fields_hits, delimiter='\t')
        fout_hits.writeheader()
        data = sorted(DictReader(open(isoform_hits_name+"_tmp"), delimiter='\t'), key=lambda row: row['Isoform'])
        for r in data:
            r['Matching_type'] = 'primary' if r['Hit'] in isoforms_info[r['Isoform']].transcripts else 'secondary'
            fout_hits.writerow(r)
    os.remove(isoform_hits_name+'_tmp')

def generate_report(saturation,report, outputClassPath, outputJuncPath):
    qc_logger.info("**** Generating SQANTI3 report.")
    cmd = f"{RSCRIPTPATH} {RSCRIPT_QC_REPORT} {outputClassPath} {outputJuncPath} {utilitiesPath} {saturation} {report}"
    logFile = f"{os.path.dirname(outputClassPath)}/logs/final_report.log"
    run_command(cmd,qc_logger,logFile,"SQANTI3 report")

def generate_tusco_report(tusco, outputClassPath, sample_gtf_file, reference_gtf_file):
    """
    Generates the TUSCO benchmarking report using the provided classification file, sample GTF, and reference GTF.
    """
    # Log to the configured handlers (StreamHandler defaults to stderr)
    qc_logger.info("Generating TUSCO report....")
    # Use species-specific TUSCO TSV (no longer depend on TUSCO GTF files)
    tusco_ref = os.path.join(utilitiesPath, "report_qc", f"tusco_{tusco}.tsv")
    if not os.path.exists(tusco_ref):
        qc_logger.error(
            f"TUSCO reference TSV not found for species '{tusco}'. Searched: {tusco_ref}"
        )
        raise FileNotFoundError("Missing TUSCO TSV reference file")
    else:
        qc_logger.info(f"Using TUSCO TSV reference: {tusco_ref}")
    # map species to genome assembly used by Gviz
    genome = {
        "human": "hg38",
        "mouse": "mm10",
    }.get(tusco, "hg38")
    # RSCRIPT_TUSCO_REPORT already contains utilitiesPath; don't prefix again
    cmd = (
        f"{RSCRIPTPATH} {RSCRIPT_TUSCO_REPORT} "
        f"{outputClassPath} {tusco_ref} {sample_gtf_file} {reference_gtf_file} {utilitiesPath} {genome}"
    )
    logFile = f"{os.path.dirname(outputClassPath)}/logs/tusco_report.log"
    # TUSCO is optional; if it fails (e.g., no matches), continue pipeline
    run_command(cmd, qc_logger, logFile, "TUSCO report", fail_ok=True)

def cleanup(outputClassPath, outputJuncPath):
    qc_logger.info("Removing temporary files.")
    try:
        os.remove(outputClassPath+"_tmp")
    except FileNotFoundError:
        pass
    try:
        os.remove(outputJuncPath+"_tmp")
    except FileNotFoundError:
        pass

    
def save_isoforms_info(isoforms_info,junctions_header, outdir, prefix):
    qc_logger.info("Saving isoforms_info object to file.")
    with open(os.path.join(outdir, f"{prefix}.isoforms_info.pkl"), 'wb') as h:
        pickle.dump(isoforms_info, h)
        pickle.dump(junctions_header, h)

def write_collapsed_GFF_with_CDS(isoforms_info, input_gff, output_gff):
    """
    Augment a collapsed GFF with CDS information
    *NEW* Also, change the "gene_id" field to use the classification result
    :param isoforms_info: dict of id -> QueryTranscript
    :param input_gff:  input GFF filename
    :param output_gff: output GFF filename
    """
    with open(output_gff, 'w') as f:
        reader = collapseGFFReader(input_gff)
        for r in reader:
            # Fix negative strand coordinates
            if r.strand == '-': 
                r.start, r.end =  r.ref_exons[-1].start, r.ref_exons[0].end # Positions get shifted by 1
                r.ref_exons.reverse() # Fix the reference exons 
            r.geneid = isoforms_info[r.seqid].geneName()  # set the gene name
            r.cds_exons = []
            if isoforms_info[r.seqid].coding == "coding": # has ORF prediction for this isoform
                s = isoforms_info[r.seqid].CDS_genomic_start  # could be 'NA'
                e = isoforms_info[r.seqid].CDS_genomic_end    # could be 'NA'
                if r.strand == '+':
                    assert s < e
                    s = s - 1 # make it 0-based
                else:
                    assert e < s
                    s, e = e, s
                    s = s - 1 # make it 0-based
                # TODO: change the loop to a binary search (reduces complexity) 
                # TODO: Include more checks into the intervals, with an equal condition
                for i,exon in enumerate(r.ref_exons):
                    if exon.end > s: break
                r.cds_exons = [Interval(s, min(e,exon.end))]
                for exon in r.ref_exons[i+1:]:
                    if exon.start > e: break
                    r.cds_exons.append(Interval(exon.start, min(e, exon.end)))
            write_collapseGFF_format(f, r)

def write_isoform_hits(isoform_hits_name,data_list):
    """
    Write isoform hits to a temporary file.
    :param isoform_hits_name: Name of the file to write the hits to.
    :param data_list: List of data to write.
    """
    with open(isoform_hits_name+'_tmp', 'a') as out_file:
        tsv_writer = writer(out_file, delimiter='\t')
        tsv_writer.writerow(data_list)


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
