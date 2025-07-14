import pickle
import shutil
import os,sys,copy,csv
import pandas as pd #type: ignore
import re

from multiprocessing import Process
from Bio import SeqIO

from src.utilities.cupcake.io.GFF import collapseGFFReader, write_collapseGFF_format

from src.config import FIELDS_CLASS
from src.qc_computations import classify_fsm, full_length_quantification, process_rts #type: ignore
from src.qc_pipeline import run
from src.helpers import get_corr_filenames, get_class_junc_filenames, get_pickle_filename, rename_novel_genes 
from src.qc_output import (
    cleanup, generate_report, write_classification_output, write_isoform_hits, write_junction_output, write_omitted_isoforms, write_collapsed_GFF_with_CDS)
from src.module_logging import qc_logger
# TODO: Create a special logging for the parallelization, to handle the individual logs of the splits into  their own files

def get_split_dir(outdir,prefix):
    split_prefix=os.path.join(os.path.abspath(outdir), prefix)
    split_directory = split_prefix+'_splits/'
    return split_directory

def natural_sort_key(s):
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', s)]

def split_input_run(args, outdir):
    SPLIT_ROOT_DIR = outdir
    if os.path.exists(SPLIT_ROOT_DIR):
        qc_logger.warning(f"{SPLIT_ROOT_DIR} directory already exists!")
    else:
        os.makedirs(SPLIT_ROOT_DIR)

    # TODO: Check here the effect of the strand when ussing collapseGFFReader
    if not args.fasta:
        try:
            recs = [r for r in collapseGFFReader(args.isoforms)]
            # Group records by gene_id
            gene_groups = {}
            for rec in recs:
                gene_id = rec.gene_id  # Assuming gene_id is an attribute of the record
                if gene_id not in gene_groups:
                    gene_groups[gene_id] = []
                gene_groups[gene_id].append(rec)
        except Exception as e:
            recs_df = pd.read_csv(args.isoforms, sep='\t', comment='#', header=None)
            for i, value in enumerate(recs_df.iloc[:, 8]):
                parts = value.split('; ')
                for part in parts:
                    if 'gene_id' in part:
                        gene_id = part.split('"')[1]
                        recs_df.at[i, 'gene_id'] = gene_id
                        break
            gene_groups = {gene_id: group.iloc[:, :-1] for gene_id, group in recs_df.groupby('gene_id')}

        n = len(gene_groups)
        if n == 0:
            qc_logger.error("The input file is not in the correct format, please check the file contains gene_id in "
                  "column 9 and try again")
            sys.exit(1)

        gene_ids = sorted(gene_groups.keys(), key=natural_sort_key)
        target_chunk_size = n // args.chunks + (1 if n % args.chunks else 0)

        split_outs = []
        current_chunk = []
        current_chunk_size = 0
        chunk_index = 0
        prev_gene_id = None

        for gene_id in gene_ids:
            if current_chunk_size >= target_chunk_size and gene_id != prev_gene_id:
                # Write the current chunk
                d = os.path.join(SPLIT_ROOT_DIR, str(chunk_index))
                os.makedirs(d, exist_ok=True)
                f = open(os.path.join(d, os.path.basename(args.isoforms) + '.split' + str(chunk_index)), 'w')
                
                for gid in current_chunk:
                    if isinstance(gene_groups[gid], pd.DataFrame):
                        gene_groups[gid].to_csv(f, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE, escapechar='\\')
                    else:
                        for rec in gene_groups[gid]:
                            write_collapseGFF_format(f, rec)
                
                f.close()
                split_outs.append((os.path.abspath(d), f.name))
                
                # Reset for the next chunk
                current_chunk = []
                current_chunk_size = 0
                chunk_index += 1

            current_chunk.append(gene_id)
            current_chunk_size += 1
            prev_gene_id = gene_id

        # Write the last chunk if it's not empty
        if current_chunk:
            d = os.path.join(SPLIT_ROOT_DIR, str(chunk_index))
            os.makedirs(d, exist_ok=True)
            f = open(os.path.join(d, os.path.basename(args.isoforms) + '.split' + str(chunk_index)), 'w')
            
            for gid in current_chunk:
                if isinstance(gene_groups[gid], pd.DataFrame):
                    gene_groups[gid].to_csv(f, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE, escapechar='\\')
                else:
                    for rec in gene_groups[gid]:
                        write_collapseGFF_format(f, rec)
            
            f.close()
            split_outs.append((os.path.abspath(d), f.name))

    else:
        # FASTA file handling remains unchanged
        recs = [r for r in SeqIO.parse(open(args.isoforms),'fasta')]
        n = len(recs)
        chunk_size = n//args.chunks + (n%args.chunks >0)
        split_outs = []
        for i in range(args.chunks):
            if i*chunk_size >= n:
                break
            d = os.path.join(SPLIT_ROOT_DIR, str(i))
            os.makedirs(d)
            f = open(os.path.join(d, os.path.basename(args.isoforms)+'.split'+str(i)), 'w')
            for j in range(i*chunk_size, min((i+1)*chunk_size, n)):
                SeqIO.write(recs[j], f, 'fasta')
            f.close()
            split_outs.append((os.path.abspath(d), f.name))

    pools = []
    for i,(d,x) in enumerate(split_outs):
        qc_logger.info(f"Launching worker on {x}....")
        args2 = copy.deepcopy(args)
        args2.isoforms = x
        args2.novel_gene_prefix = str(i)
        args2.dir = d
        args2.report = 'skip'
        p = Process(target=run, args=(args2,))
        p.start()
        pools.append(p)

    for p in pools:
        p.join()
    return [d for (d,x) in split_outs]

def combine_split_runs(args, split_dirs):
    """
    Combine .faa, .fasta, .gtf, .classification.txt, .junctions.txt
    Then write out the PDF report
    """
    corrGTF, _, corrFASTA, corrORF , corrCDS_GTF_GFF = get_corr_filenames(args.dir, args.output)
    outputClassPath, outputJuncPath = get_class_junc_filenames(args.dir,args.output)

    if not args.skipORF:
        f_faa = open(corrORF, 'w')
    f_fasta = open(corrFASTA, 'w')
    f_gtf = open(corrGTF, 'w')
    f_junc_temp = open(outputJuncPath+"_tmp", "w")
    isoforms_info = {}
    headers = []
    for i,split_d in enumerate(split_dirs):
        _gtf, _, _fasta, _orf , _ = get_corr_filenames(split_d,args.output)
        _, _junc = get_class_junc_filenames(split_d,args.output)
        _info = get_pickle_filename(split_d,args.output)
        if not args.skipORF:
            with open(_orf) as h: f_faa.write(h.read())
        with open(_gtf) as h: f_gtf.write(h.read())
        with open(_fasta) as h: f_fasta.write(h.read())
        with open(f"{_junc}_tmp") as h:
            if i == 0:
                f_junc_temp.write(h.readline())
            else:
                h.readline()
            f_junc_temp.write(h.read())
        # Retrieving the isoform object and the junctions header
        with open(_info, 'rb') as h:
            isoforms_info.update(pickle.load(h))
            headers.append(pickle.load(h))
        shutil.move(os.path.join(split_d,"TD2"),os.path.join(args.dir,"TD2",f"TD2_{i}"))
    f_fasta.close()
    f_gtf.close()
    f_junc_temp.close()
    
    # Fix novel genes and classify FSM
    isoforms_info = rename_novel_genes(isoforms_info, args.novel_gene_prefix)
    isoforms_info = classify_fsm(isoforms_info)
    ## FL count file
    fields_class_cur = FIELDS_CLASS
    if args.fl_count:
        isoforms_info, fields_class_cur = full_length_quantification(args.fl_count, isoforms_info, FIELDS_CLASS)
    isoforms_info,RTS_info = process_rts(isoforms_info,outputJuncPath,args.refFasta)

    fields_junc_cur = headers[0]
    if not args.skipORF:
        write_collapsed_GFF_with_CDS(isoforms_info, corrGTF, corrCDS_GTF_GFF)
    write_classification_output(isoforms_info, outputClassPath, FIELDS_CLASS)
    write_junction_output(outputJuncPath, RTS_info, fields_junc_cur)
    #write omitted isoforms if requested minimum reference length is more than 0
    isoforms_info = write_omitted_isoforms(isoforms_info, args.dir, args.output, 
                                            args.min_ref_len, args.is_fusion, fields_class_cur)
    cleanup(outputClassPath, outputJuncPath)
    #isoform hits to file if requested
    if args.isoform_hits:
        write_isoform_hits(args.dir, args.output, isoforms_info)

    if not args.skipORF:
        f_faa.close()

    if args.report != 'skip':
        generate_report(args.saturation,args.report, outputClassPath, outputJuncPath)
