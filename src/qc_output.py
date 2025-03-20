from csv import DictReader, DictWriter
import os
import pickle
import sys

from bx.intervals import Interval #type: ignore

from src.utilities.cupcake.io.GFF import collapseGFFReader, write_collapseGFF_format

from .commands import (
    RSCRIPTPATH, RSCRIPT_REPORT, run_command,
    utilitiesPath
)
from src.helpers import get_isoform_hits_name, get_omitted_name


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
     
    # iso_keys = sorted(isoforms_info.keys(), key=lambda x: (isoforms_info[x].chrom, isoforms_info[x].id))
    with open(outputClassPath, 'w') as h:
        fout_class = DictWriter(h, fieldnames=fields_class_cur, delimiter='\t')
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
    print("**** Generating SQANTI3 report....", file=sys.stderr)
    cmd = f"{RSCRIPTPATH} {utilitiesPath}/{RSCRIPT_REPORT} {outputClassPath} {outputJuncPath} {utilitiesPath} {saturation} {report}"
    run_command(cmd,"SQANTI3 report")

def cleanup(outputClassPath, outputJuncPath):
    print("Removing temporary files....", file=sys.stderr)
    try:
        os.remove(outputClassPath+"_tmp")
    except FileNotFoundError:
        pass
    try:
        os.remove(outputJuncPath+"_tmp")
    except FileNotFoundError:
        pass
    
def save_isoforms_info(isoforms_info,junctions_header, outdir, prefix):
    print("Saving isoforms_info object to file....", file=sys.stderr)
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
            print(isoforms_info[r.seqid].genes)
            r.geneid = isoforms_info[r.seqid].geneName()  # set the gene name
            s = isoforms_info[r.seqid].CDS_genomic_start  # could be 'NA'
            e = isoforms_info[r.seqid].CDS_genomic_end    # could be 'NA'
            r.cds_exons = []
            if s!='NA' and e!='NA': # has ORF prediction for this isoform
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