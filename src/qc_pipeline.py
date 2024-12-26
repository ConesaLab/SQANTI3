import os, timeit, sys
import subprocess

from csv import DictReader, DictWriter
from Bio import SeqIO

from src.utilities.indels_annot import calc_indels_from_sam

from .helpers import (
    get_corr_filenames, get_class_junc_filenames, get_isoform_hits_name,
    get_omitted_name,sequence_correction, predictORF, write_collapsed_GFF_with_CDS
    )
from .parsers import (
    reference_parser, isoforms_parser
)
from .classification import isoformClassification, preprocess_isoform_data
from .config import FIELDS_CLASS 
from .commands import (
    RSCRIPTPATH, RSCRIPT_REPORT, ISOANNOT_PROG,
    utilitiesPath, GTF_to_genePred
)
from src.qc_computations import (
    classify_fsm, isoform_expression_info, isoforms_junctions,
    process_rts_swiching, ratio_TSS_dict_reading,
    full_length_quantification
)

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
    
    orfDict = predictORF(args.dir, args.skipORF, args.orf_input, 
                         corrFASTA, corrORF)

    ## parse reference id (GTF) to dicts
    refs_1exon_by_chr, refs_exons_by_chr, \
        junctions_by_chr, junctions_by_gene, start_ends_by_gene = \
            reference_parser(args.annotation,args.dir,args.output, list(genome_dict.keys()),
                             args.genename, args.isoAnnotLite)

    ## parse query isoforms
    corrgenPred = GTF_to_genePred(corrGTF)
    isoforms_by_chr = isoforms_parser(corrgenPred)

    ## Run indel computation if sam exists
    # indelsJunc: dict of pbid --> list of junctions near indel (in Interval format)
    # indelsTotal: dict of pbid --> total indels count
    if os.path.exists(corrSAM):
        (indelsJunc, indelsTotal) = calc_indels_from_sam(corrSAM)
    else:
        indelsJunc = None
        indelsTotal = None
    
    fusion_components, isoform_hits_name, SJcovNames, SJcovInfo, \
    fields_junc_cur,ratio_TSS_dict, cage_peak_obj, polya_peak_obj, \
    polyA_motif_list, phyloP_reader = preprocess_isoform_data(args, corrGTF)
    
    # isoform classification + intra-priming + id and junction characterization
    isoforms_info, ratio_TSS_dict = isoformClassification(args.sites,args.window,args.novel_gene_prefix,args.is_fusion,
                                                          isoforms_by_chr, refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr, 
                          junctions_by_gene, start_ends_by_gene, genome_dict, indelsJunc, orfDict,
                          outputClassPath, outputJuncPath,fusion_components, isoform_hits_name,SJcovNames, SJcovInfo,
                          fields_junc_cur,ratio_TSS_dict, cage_peak_obj, polya_peak_obj, polyA_motif_list, phyloP_reader)

    print("Number of classified isoforms: {0}".format(len(isoforms_info)), file=sys.stdout)

    write_collapsed_GFF_with_CDS(isoforms_info, corrGTF, corrGTF+'.cds.gff')

    ## RT-switching computation
    print("**** RT-switching computation....", file=sys.stderr)

    # RTS_info: dict of (pbid) -> list of RT junction. if RTS_info[pbid] == [], means all junctions are non-RT.
    
    isoforms_info, RTS_info = process_rts_swiching(isoforms_info,outputJuncPath,
                                                    args.genome,genome_dict)
    print(f"After RTS classificaion: {len(isoforms_info)}")

    ## FSM classification
    isoforms_info = classify_fsm(isoforms_info)
    print(f"After classify fsm: {len(isoforms_info)}")

    fields_class_cur = FIELDS_CLASS
    ## FL count file
    if args.fl_count:
        isoforms_info, fields_class_cur = full_length_quantification(args.fl_count, isoforms_info, FIELDS_CLASS)
    else:
        print("Full-length read abundance files not provided.", file=sys.stderr)
    
    ## TSS ratio dict reading
    if ratio_TSS_dict is not None:
        isoforms_info = ratio_TSS_dict_reading(isoforms_info, ratio_TSS_dict)

    ## Isoform expression information
    isoforms_info = isoform_expression_info(isoforms_info,args.expression,args.short_reads,
                                           args.dir,corrFASTA,args.cpus)

    if indelsTotal is not None:
        for iso in isoforms_info:
            if iso in indelsTotal:
                isoforms_info[iso].nIndels = indelsTotal[iso]
            else:
                isoforms_info[iso].nIndels = 0


    ## Read junction files and create attributes per id
    reader = DictReader(open(outputJuncPath+"_tmp"), delimiter='\t')
    fields_junc_cur = reader.fieldnames
    isoforms_info = isoforms_junctions(isoforms_info, reader)

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
    os.remove(outputClassPath+"_tmp")
    os.remove(outputJuncPath+"_tmp")

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


