import os, timeit, sys
import subprocess

from csv import DictReader
from Bio import SeqIO  # type: ignore

from src.utilities.indels_annot import calc_indels_from_sam

from src.qc_output import (
    cleanup, generate_report, generate_tusco_report, save_isoforms_info,
    write_classification_output, write_isoform_hits, write_junction_output,
    write_omitted_isoforms, write_collapsed_GFF_with_CDS
)
from src.helpers import (
    get_corr_filenames, get_class_junc_filenames, rename_novel_genes, 
    sequence_correction, predictORF
    )
from src.parsers import (
    get_fusion_component, reference_parser, isoforms_parser
)
from src.config import FIELDS_CLASS 
from src.commands import (
    ISOANNOT_PROG, GTF_to_genePred
)
from src.qc_computations import (
    classify_fsm, isoform_expression_info, isoforms_junctions,
    process_rts, ratio_TSS_dict_reading,
    full_length_quantification
)
from src.classification_preprocessing import (
    initialize_isoform_hits, read_CAGE_peaks, read_polyA_peaks,
    read_polyA_motifs, read_phyloP_bed, SJ_coverage, TSS_ratio_calculation
)
from src.classification_main import isoform_classification_pipeline
from src.module_logging import qc_logger

def run(args):
    global isoform_hits_name

    # Get filenames
    corrGTF, corrSAM, corrFASTA, corrORF , corrCDS_GTF_GFF = get_corr_filenames(args.dir, args.output) # Fix the path to the corrected CDS file
    badstrandGTF = args.dir + "/unknown_strand.gtf"
    outputClassPath, outputJuncPath = get_class_junc_filenames(args.dir,args.output)

    start3 = timeit.default_timer()

    qc_logger.info("Parsing provided files")
    qc_logger.info(f"Reading genome fasta {args.refFasta}")
    # NOTE: can't use LazyFastaReader because inefficient. Bring the whole genome in!
    genome_dict = dict((r.name, r) for r in SeqIO.parse(open(args.refFasta), 'fasta'))
    ## correction of sequences and ORF prediction (if gtf provided instead of fasta file, correction of sequences will be skipped)
    sequence_correction(
    args.dir, args.output, args.cpus, args.chunks, args.fasta,
    genome_dict, badstrandGTF, args.refFasta, args.isoforms, args.aligner_choice,
    gmap_index=args.gmap_index, annotation=args.refGTF)
    
    orfDict = predictORF(args.dir, args.skipORF, args.orf_input, 
                         corrFASTA, corrORF,args.cpus)

    ## parse reference id (GTF) to dicts
    refs_1exon_by_chr, refs_exons_by_chr, \
        junctions_by_chr, junctions_by_gene, start_ends_by_gene = \
            reference_parser(args.refGTF,args.dir,args.output, list(genome_dict.keys()),
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
    
    # Preprocess isoform data
    fusion_components =  get_fusion_component(args.isoforms) if args.is_fusion else {}
    isoform_hits_name = initialize_isoform_hits(args.dir, args.output,args.isoform_hits)
    star_out, star_index, SJcovNames,\
        SJcovInfo, fields_junc_cur = SJ_coverage(args.short_reads, args.coverage, 
                                                 args.refFasta, args.dir, args.cpus)
    ratio_TSS_dict = TSS_ratio_calculation(args.SR_bam, args.short_reads, star_out,
                                           star_index, corrGTF, args.ratio_TSS_metric)
    cage_peak_obj = read_CAGE_peaks(args.CAGE_peak)
    polya_peak_obj = read_polyA_peaks(args.polyA_peak)
    polyA_motif_list = read_polyA_motifs(args.polyA_motif_list)
    phyloP_reader = read_phyloP_bed(args.phyloP_bed)

    # isoform classification + intra-priming + id and junction characterization
    isoforms_info, ratio_TSS_dict = isoform_classification_pipeline(
        args.sites, args.window, args.is_fusion,
        isoforms_by_chr, refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr,
        junctions_by_gene, start_ends_by_gene, genome_dict, indelsJunc, orfDict,
        outputClassPath, outputJuncPath, fusion_components, isoform_hits_name,
        SJcovNames, SJcovInfo, fields_junc_cur, ratio_TSS_dict, cage_peak_obj,
        polya_peak_obj, polyA_motif_list, phyloP_reader)


    qc_logger.info(f"Number of classified isoforms: {len(isoforms_info)}")
    
    # This steps are avoided in the parallel implementation
    # They are run instead after combining the chunks
    if args.chunks == 1:
        ## FSM classification
        isoforms_info = rename_novel_genes(isoforms_info, args.novel_gene_prefix)
        isoforms_info = classify_fsm(isoforms_info)
        qc_logger.info(f"After classify fsm: {len(isoforms_info)}")
        
        ## FL count file
        fields_class_cur = FIELDS_CLASS
        if args.fl_count:
            isoforms_info, fields_class_cur = full_length_quantification(args.fl_count, isoforms_info, FIELDS_CLASS)
        else:
            qc_logger.info("Full-length read abundance files not provided.")

        ## RT-switching computation
        qc_logger.info("RT-switching computation")
        isoforms_info, RTS_info = process_rts(isoforms_info,outputJuncPath,
                                                        args.refFasta,genome_dict)
        qc_logger.info(f"After RTS classificaion: {len(isoforms_info)}")

    
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

    #### qc_logger.infoing output file:
    qc_logger.info("Writing output files")

    if args.chunks != 1:
        save_isoforms_info(isoforms_info, fields_junc_cur, args.dir, args.output)
        
    else:
        # Write corrected GTF with CDS
        if not args.skipORF:
            write_collapsed_GFF_with_CDS(isoforms_info, corrGTF, corrCDS_GTF_GFF)

        # Write final classification
        write_classification_output(isoforms_info, outputClassPath, fields_class_cur)

        # Now that RTS info is obtained, we can write the final junctions.txt
        write_junction_output(outputJuncPath, RTS_info, fields_junc_cur)

        #write omitted isoforms if requested minimum reference length is more than 0
        isoforms_info = write_omitted_isoforms(isoforms_info, args.dir, args.output, 
                                            args.min_ref_len, args.is_fusion, fields_class_cur)
        
        #isoform hits to file if requested
        if args.isoform_hits:
            write_isoform_hits(args.dir, args.output, isoforms_info)
        
        ## Generating report
        if args.report != 'skip':
            # Run TUSCO benchmarking report if requested
            if hasattr(args, 'tusco') and args.tusco:
                # Pass both sample GTF (corrected) and reference GTF for IGV-like plots
                generate_tusco_report(args.tusco, outputClassPath, corrGTF, args.refGTF)
            # Main SQANTI3 report
            generate_report(args.saturation, args.report, outputClassPath, outputJuncPath)

        cleanup(outputClassPath, outputJuncPath)

    stop3 = timeit.default_timer()
    qc_logger.info(f"SQANTI3 complete in {stop3 - start3} sec.")


### IsoAnnot Lite implementation
# ISOANNOT_PROG =  os.path.join(utilitiesPath, "IsoAnnotLite_SQ3.py")

# TODO: Change this so it is part of SQANTI, called as a function and not as a external script
def run_isoAnnotLite(correctedGTF, outClassFile, outJuncFile, outName, gff3_opt):
    if gff3_opt:
        iso_out = os.path.join(os.path.dirname(correctedGTF), outName)
        isoAnnot_sum = iso_out + ".isoAnnotLite_stats.txt"
        ISOANNOT_CMD = f"python3 {ISOANNOT_PROG} {correctedGTF} {outClassFile} {outJuncFile} -gff3 {gff3_opt} -o {iso_out} -novel -stdout {isoAnnot_sum}"
    else:
        iso_out = os.path.join(os.path.dirname(correctedGTF), outName)
        ISOANNOT_CMD = f"python3 {ISOANNOT_PROG} {correctedGTF} {outClassFile} {outJuncFile} -o {iso_out} -novel"
    if subprocess.check_call(ISOANNOT_CMD, shell=True)!=0:
        qc_logger.error(f"Command failed: {ISOANNOT_CMD}")
        sys.exit(1)
