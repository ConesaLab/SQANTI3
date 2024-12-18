import os, timeit, sys
import subprocess

from csv import DictReader, DictWriter
from collections import defaultdict
from Bio import SeqIO

import pandas as pd
from src.utilities.rt_switching import rts
from src.utilities.indels_annot import calc_indels_from_sam
from src.utilities.short_reads import kallisto

from .helpers import (
    get_corr_filenames, get_class_junc_filenames, get_isoform_hits_name,
    get_omitted_name,sequence_correction, predictORF, write_collapsed_GFF_with_CDS
    )
from .parsers import reference_parser, isoforms_parser, FLcount_parser, expression_parser
from .classification import isoformClassification
from .config import FIELDS_CLASS 
from .commands import RSCRIPTPATH, RSCRIPT_REPORT, ISOANNOT_PROG, utilitiesPath, short_reads_mapping
from .utils import pstdev

def run(args):

    global isoform_hits_name

    corrGTF, corrSAM, corrFASTA, corrORF , _ = get_corr_filenames(args.dir, args.output)
    badstrandGTF = args.dir + "/unknown_strand.gtf"
    outputClassPath, outputJuncPath = get_class_junc_filenames(args.dir,args.output)

    start3 = timeit.default_timer()

    print("**** Parsing provided files....", file=sys.stdout)
    print("Reading genome fasta {0}....".format(args.genome), file=sys.stdout)
    # NOTE: can't use LazyFastaReader because inefficient. Bring the whole genome in!
    genome_dict = dict((r.name, r) for r in SeqIO.parse(open(args.genome), 'fasta'))

    ## correction of sequences and ORF prediction (if gtf provided instead of fasta file, correction of sequences will be skipped)
    sequence_correction(
    args.dir, args.output, args.cpus, args.chunks, args.fasta,
    genome_dict, badstrandGTF, args.genome, args.isoforms, args.aligner_choice,
    gmap_index=args.gmap_index, sense=args.sense, annotation=args.annotation)
    
    if not os.path.exists(corrFASTA):
        print("ERROR: corrected FASTA file {0} does not exist! Abort!".format(corrFASTA), file=sys.stderr)
        sys.exit(1)
    else:
        print("Corrected FASTA file written to {0}.".format(corrFASTA), file=sys.stdout)
    orfDict = predictORF(args, corrFASTA, corrORF)

    ## parse reference id (GTF) to dicts
    refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr, junctions_by_gene, start_ends_by_gene = reference_parser(args.annotation,args.dir,args.output,
                                                                                                                     args.genename,                                                                                                                     args.isoAnnotLite, 
                                                                                                                     list(genome_dict.keys()))

    ## parse query isoforms
    isoforms_by_chr = isoforms_parser(corrGTF)

    ## Run indel computation if sam exists
    # indelsJunc: dict of pbid --> list of junctions near indel (in Interval format)
    # indelsTotal: dict of pbid --> total indels count
    if os.path.exists(corrSAM):
        (indelsJunc, indelsTotal) = calc_indels_from_sam(corrSAM)
    else:
        indelsJunc = None
        indelsTotal = None

    ## Short-read mapping
    star_out, star_index, SJcovNames, SJcovInfo = short_reads_mapping(args)

    # isoform classification + intra-priming + id and junction characterization
    isoforms_info, ratio_TSS_dict = isoformClassification(args, isoforms_by_chr, refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr, junctions_by_gene, start_ends_by_gene, genome_dict, indelsJunc, orfDict, 
                                                          corrGTF, star_out, star_index, SJcovNames, SJcovInfo,outputClassPath, outputJuncPath )

    print("Number of classified isoforms: {0}".format(len(isoforms_info)), file=sys.stdout)

    write_collapsed_GFF_with_CDS(isoforms_info, corrGTF, corrGTF+'.cds.gff')

    ## RT-switching computation
    print("**** RT-switching computation....", file=sys.stderr)

    # RTS_info: dict of (pbid) -> list of RT junction. if RTS_info[pbid] == [], means all junctions are non-RT.
    RTS_info = rts([outputJuncPath+"_tmp", args.genome, "-a"], genome_dict)
    for pbid in isoforms_info:
        if pbid in RTS_info and len(RTS_info[pbid]) > 0:
            isoforms_info[pbid].RT_switching = "TRUE"
        else:
            isoforms_info[pbid].RT_switching = "FALSE"


    ## FSM classification
    geneFSM_dict = defaultdict(lambda: [])
    for iso in isoforms_info:
        gene = isoforms_info[iso].geneName()  # if multi-gene, returns "geneA_geneB_geneC..."
        geneFSM_dict[gene].append(isoforms_info[iso].str_class)

    fields_class_cur = FIELDS_CLASS
    ## FL count file
    if args.fl_count:
        if not os.path.exists(args.fl_count):
            print("FL count file {0} does not exist!".format(args.fl_count), file=sys.stderr)
            sys.exit(1)
        print("**** Reading Full-length read abundance files...", file=sys.stderr)
        fl_samples, fl_count_dict = FLcount_parser(args.fl_count)
        for pbid in fl_count_dict:
            if pbid not in isoforms_info:
                print("WARNING: {0} found in FL count file but not in input fasta.".format(pbid), file=sys.stderr)
        if len(fl_samples) == 1: # single sample from PacBio
            print("Single-sample PacBio FL count format detected.", file=sys.stderr)
            for iso in isoforms_info:
                if iso in fl_count_dict:
                    isoforms_info[iso].FL = fl_count_dict[iso]
                else:
                    print("WARNING: {0} not found in FL count file. Assign count as 0.".format(iso), file=sys.stderr)
                    isoforms_info[iso].FL = 0
        else: # multi-sample
            print("Multi-sample PacBio FL count format detected.", file=sys.stderr)
            fields_class_cur = FIELDS_CLASS + ["FL."+s for s in fl_samples]
            for iso in isoforms_info:
                if iso in fl_count_dict:
                    isoforms_info[iso].FL_dict = fl_count_dict[iso]
                else:
                    print("WARNING: {0} not found in FL count file. Assign count as 0.".format(iso), file=sys.stderr)
                    isoforms_info[iso].FL_dict = defaultdict(lambda: 0)
    else:
        print("Full-length read abundance files not provided.", file=sys.stderr)

    ## TSS ratio dict reading
    if ratio_TSS_dict is not None:
        print('**** Adding TSS ratio data... ****')
        for iso in ratio_TSS_dict:
            if iso not in isoforms_info:
                print("WARNING: {0} found in ratio TSS file but not in input FASTA/GTF".format(iso), file=sys.stderr)
        for iso in isoforms_info:
            if iso in ratio_TSS_dict:
                if str(ratio_TSS_dict[iso]['return_ratio']) == 'nan':
                    isoforms_info[iso].ratio_TSS = 'NA'
                else:
                    isoforms_info[iso].ratio_TSS = ratio_TSS_dict[iso]['return_ratio']
            else:
                print("WARNING: {0} not found in ratio TSS file. Assign count as 1.".format(iso), file=sys.stderr)
                isoforms_info[iso].ratio_TSS = 1

    ## Isoform expression information
    if args.expression:
        print("**** Reading Isoform Expression Information.", file=sys.stderr)
        exp_dict = expression_parser(args.expression)
        gene_exp_dict = {}
        for iso in isoforms_info:
            if iso not in exp_dict:
                exp_dict[iso] = 0
                print("WARNING: isoform {0} not found in expression matrix. Assigning TPM of 0.".format(iso), file=sys.stderr)
            gene = isoforms_info[iso].geneName()
            if gene not in gene_exp_dict:
                gene_exp_dict[gene] = exp_dict[iso]
            else:
                gene_exp_dict[gene] = gene_exp_dict[gene]+exp_dict[iso]
    else:
        if args.short_reads is not None:
            print("**** Running Kallisto to calculate isoform expressions. ")
            expression_files = kallisto(corrFASTA, args.short_reads, args.dir, args.cpus)
            exp_dict = expression_parser(expression_files)
            gene_exp_dict = {}
            for iso in isoforms_info:
                if iso not in exp_dict:
                    exp_dict[iso] = 0
                    print("WARNING: isoform {0} not found in expression matrix. Assigning TPM of 0.".format(iso), file=sys.stderr)
                gene = isoforms_info[iso].geneName()
                if gene not in gene_exp_dict:
                    gene_exp_dict[gene] = exp_dict[iso]
                else:
                    gene_exp_dict[gene] = gene_exp_dict[gene]+exp_dict[iso]
        else:
            exp_dict = None
            gene_exp_dict = None
            print("Isoforms expression files not provided.", file=sys.stderr)


    ## Adding indel, FSM class and expression information
    for iso in isoforms_info:
        gene = isoforms_info[iso].geneName()
        if exp_dict is not None and gene_exp_dict is not None:
            isoforms_info[iso].geneExp = gene_exp_dict[gene]
            isoforms_info[iso].isoExp  = exp_dict[iso]
        if len(geneFSM_dict[gene])==1:
            isoforms_info[iso].FSM_class = "A"
        elif "full-splice_match" in geneFSM_dict[gene]:
            isoforms_info[iso].FSM_class = "C"
        else:
            isoforms_info[iso].FSM_class = "B"

    if indelsTotal is not None:
        for iso in isoforms_info:
            if iso in indelsTotal:
                isoforms_info[iso].nIndels = indelsTotal[iso]
            else:
                isoforms_info[iso].nIndels = 0


    ## Read junction files and create attributes per id
    # Read the junction information to fill in several remaining unfilled fields in classification
    # (1) "canonical": is "canonical" if all junctions are canonical, otherwise "non_canonical"
    # (2) "bite": is TRUE if any of the junction "bite_junction" field is TRUE

    reader = DictReader(open(outputJuncPath+"_tmp"), delimiter='\t')
    fields_junc_cur = reader.fieldnames

    sj_covs_by_isoform = defaultdict(lambda: [])  # pbid --> list of total_cov for each junction so we can calculate SD later
    for r in reader:
        # only need to do assignment if:
        # (1) the .canonical field is still "NA"
        # (2) the junction is non-canonical
        assert r['canonical'] in ('canonical', 'non_canonical')
        if (isoforms_info[r['isoform']].canonical == 'NA') or \
            (r['canonical'] == 'non_canonical'):
            isoforms_info[r['isoform']].canonical = r['canonical']

        if (isoforms_info[r['isoform']].bite == 'NA') or (r['bite_junction'] == 'TRUE'):
            isoforms_info[r['isoform']].bite = r['bite_junction']

        if r['indel_near_junct'] == 'TRUE':
            if isoforms_info[r['isoform']].nIndelsJunc == 'NA':
                isoforms_info[r['isoform']].nIndelsJunc = 0
            isoforms_info[r['isoform']].nIndelsJunc += 1

        # min_cov: min( total_cov[j] for each junction j in this isoform )
        # min_cov_pos: the junction [j] that attributed to argmin(total_cov[j])
        # min_sample_cov: min( sample_cov[j] for each junction in this isoform )
        # sd_cov: sd( total_cov[j] for each junction j in this isoform )
        if r['sample_with_cov'] != 'NA':
            sample_with_cov = int(r['sample_with_cov'])
            if (isoforms_info[r['isoform']].min_samp_cov == 'NA') or (isoforms_info[r['isoform']].min_samp_cov > sample_with_cov):
                isoforms_info[r['isoform']].min_samp_cov = sample_with_cov

        if r['total_coverage_unique'] != 'NA':
            total_cov = int(r['total_coverage_unique'])
            sj_covs_by_isoform[r['isoform']].append(total_cov)
            if (isoforms_info[r['isoform']].min_cov == 'NA') or (isoforms_info[r['isoform']].min_cov > total_cov):
                isoforms_info[r['isoform']].min_cov = total_cov
                isoforms_info[r['isoform']].min_cov_pos = r['junction_number']


    for pbid, covs in sj_covs_by_isoform.items():
        isoforms_info[pbid].sd = pstdev(covs)

    #### Printing output file:
    print("**** Writing output files....", file=sys.stderr)

    #write omitted isoforms if requested minimum reference length is more than 0
    if args.min_ref_len > 0 and not args.is_fusion:
        omitted_name = get_omitted_name(args.dir, args.output)
        omitted_iso = {}
        for key in isoforms_info:
            if not isoforms_info[key].refLen == 'NA':
                if int(isoforms_info[key].refLen) <= int(args.min_ref_len):
                    omitted_iso[key] = isoforms_info[key]
        for key in omitted_iso:
            del isoforms_info[key]
        omitted_keys = list(omitted_iso.keys())
        # TODO: check if this print is necessary
        print(type(omitted_keys))
        omitted_keys.sort(key=lambda x: (omitted_iso[x].chrom,omitted_iso[x].id))
        print('Type omitted keys ', type(omitted_keys))
        with open(omitted_name, 'w') as h:
            fout_omitted = DictWriter(h, fieldnames=fields_class_cur, delimiter='\t')
            fout_omitted.writeheader()
            for key in omitted_keys:
                fout_omitted.writerow(omitted_iso[key].as_dict())
    # sort isoform keys
    iso_keys = list(isoforms_info.keys())
    iso_keys.sort(key=lambda x: (isoforms_info[x].chrom,isoforms_info[x].id))
    with open(outputClassPath, 'w') as h:
        fout_class = DictWriter(h, fieldnames=fields_class_cur, delimiter='\t')
        fout_class.writeheader()
        for iso_key in iso_keys:
            fout_class.writerow(isoforms_info[iso_key].as_dict())

    # Now that RTS info is obtained, we can write the final junctions.txt
    with open(outputJuncPath, 'w') as h:
        fout_junc = DictWriter(h, fieldnames=fields_junc_cur, delimiter='\t')
        fout_junc.writeheader()
        for r in DictReader(open(outputJuncPath+"_tmp"), delimiter='\t'):
            if r['isoform'] in RTS_info:
                if r['junction_number'] in RTS_info[r['isoform']]:
                    r['RTS_junction'] = 'TRUE'
                else:
                    r['RTS_junction'] = 'FALSE'
            fout_junc.writerow(r)
    #isoform hits to file if requested
    if args.isoform_hits:
        fields_hits =['Isoform', 'Isoform_length', 'Isoform_exon_number', 'Hit', 'Hit_length',
                      'Hit_exon_number', 'Match', 'Diff_to_TSS', 'Diff_to_TTS', 'Matching_type']
        isoform_hits_name = get_isoform_hits_name(args.dir, args.output)
        with open(isoform_hits_name,'w') as h:
            fout_hits = DictWriter(h, fieldnames=fields_hits, delimiter='\t')
            fout_hits.writeheader()
            data = DictReader(open(isoform_hits_name+"_tmp"), delimiter='\t')
            data = sorted(data, key=lambda row:(row['Isoform']))
            for r in data:
                if r['Hit'] in isoforms_info[r['Isoform']].transcripts:
                    r['Matching_type'] = 'primary'
                else:
                    r['Matching_type'] = 'secondary'
                fout_hits.writerow(r)
        os.remove(isoform_hits_name+'_tmp')
    ## Generating report
    if args.report != 'skip':
        print("**** Generating SQANTI3 report....", file=sys.stderr)
        cmd = RSCRIPTPATH + " {d}/{f} {c} {j} {p} {d} {a} {b}".format(d=utilitiesPath, f=RSCRIPT_REPORT, c=outputClassPath, j=outputJuncPath, p=args.doc, a=args.saturation, b=args.report)
        if subprocess.check_call(cmd, shell=True)!=0:
            print("ERROR running command: {0}".format(cmd), file=sys.stderr)
            sys.exit(1)
    stop3 = timeit.default_timer()

    print("Removing temporary files....", file=sys.stderr)
    #os.remove(outputClassPath+"_tmp")
    #os.remove(outputJuncPath+"_tmp")

    print("SQANTI3 complete in {0} sec.".format(stop3 - start3), file=sys.stderr)


### IsoAnnot Lite implementation
# ISOANNOT_PROG =  os.path.join(utilitiesPath, "IsoAnnotLite_SQ3.py")

def run_isoAnnotLite(correctedGTF, outClassFile, outJuncFile, outName, gff3_opt):
    if gff3_opt:
        iso_out = os.path.join(os.path.dirname(correctedGTF), outName)
        isoAnnot_sum = iso_out + ".isoAnnotLite_stats.txt"
        ISOANNOT_CMD = "python3 "+ ISOANNOT_PROG + " {g} {c} {j} -gff3 {t} -o {o} -novel -stdout {i}".format(g=correctedGTF , c=outClassFile, j=outJuncFile, t=gff3_opt, o=iso_out, i=isoAnnot_sum)
    else:
        iso_out = os.path.join(os.path.dirname(correctedGTF), outName)
        ISOANNOT_CMD = "python3 "+ ISOANNOT_PROG + " {g} {c} {j} -o {o} -novel".format(g=correctedGTF , c=outClassFile, j=outJuncFile, o=iso_out)
    if subprocess.check_call(ISOANNOT_CMD, shell=True)!=0:
        print("ERROR running command: {0}".format(ISOANNOT_CMD), file=sys.stderr)
        sys.exit(1)


