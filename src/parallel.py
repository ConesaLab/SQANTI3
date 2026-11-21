import pickle
import shutil
import os,copy
import re
import tempfile
import time
import psutil

from multiprocessing import Process, get_start_method
from multiprocessing.connection import wait
from Bio import SeqIO


from src.qc_computations import classify_fsm, full_length_quantification, process_rts, isoform_expression_info #type: ignore
from src.qc_pipeline import run
from src.parsers import reference_parser
from src.helpers import (
    align_isoforms, get_corr_filenames, get_class_junc_filenames, get_isoform_hits_name, get_pickle_filename, rename_novel_genes, 
    get_supplementary_name, get_unmapped_name, is_fastq_file, read_fasta_fastq, uncompressed_size, warn_supplementary, warn_unmapped)
from src.qc_output import (
    cleanup, generate_report, generate_tusco_report, write_classification_output, write_isoform_hits, write_junction_output, write_omitted_isoforms, write_collapsed_GFF_with_CDS)
from src.module_logging import qc_logger

# TODO: Create a special logging for the parallelization, to handle the individual logs of the splits into  their own files

def get_split_dir(outdir,prefix):
    split_prefix=os.path.join(os.path.abspath(outdir), prefix)
    split_directory = split_prefix+'_splits/'
    return split_directory

def natural_sort_key(s):
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', s)]

def split_fasta_input(isoforms, split_root_dir, chunks):
    """
    Split a FASTA or FASTQ file into at most `chunks` FASTA files, one per directory in split_root_dir.

    :return: list of (chunk directory, chunk FASTA file)
    """
    # FASTA or FASTQ input, streamed in one pass (read files can be very large). Sequences are written
    # in input order, and a new chunk starts once ~file_size/chunks bytes of input have been read,
    # as in the GTF branch. Chunks are written as FASTA: FASTQ qualities are not used downstream.
    # Gzipped input uses its uncompressed size, so it is split exactly like the decompressed file.
    file_size = uncompressed_size(isoforms)
    target_chunk_bytes = max(1, file_size // chunks + (1 if file_size % chunks else 0))
    is_fastq = is_fastq_file(isoforms)

    split_outs = []
    current_f = None
    current_chunk_bytes = 0
    try:
        for title, seq in read_fasta_fastq(isoforms):
            if current_f is None or (current_chunk_bytes >= target_chunk_bytes and len(split_outs) < chunks):
                if current_f is not None:
                    current_f.close()
                d = os.path.join(split_root_dir, str(len(split_outs)))
                os.makedirs(d, exist_ok=True)
                current_f = open(os.path.join(d, os.path.basename(isoforms) + '.split' + str(len(split_outs))), 'w')
                split_outs.append((os.path.abspath(d), current_f.name))
                current_chunk_bytes = 0
            current_f.write(f">{title}\n{seq}\n")
            # Size of the record in the input file: header, sequence (+ '+' line and qualities in FASTQ)
            current_chunk_bytes += len(title) + 2 + (2 * len(seq) + 4 if is_fastq else len(seq) + 1)
    finally:
        if current_f is not None and not current_f.closed:
            current_f.close()

    if not split_outs:
        qc_logger.error(f"Input file '{isoforms}' contains no sequences.")
        raise ValueError(f"Input file '{isoforms}' contains no sequences.")
    return split_outs

# Root, limit, usage and reclaimable-cache names of each cgroup memory layout
CGROUP_MEMORY_LAYOUTS = {
    'memory': ('/sys/fs/cgroup/memory', 'memory.limit_in_bytes',
               'memory.usage_in_bytes', 'total_inactive_file'),  # cgroup v1
    '': ('/sys/fs/cgroup', 'memory.max', 'memory.current',
         'inactive_file'),  # cgroup v2
}

def read_cgroup_value(path):
    """Return the integer in a cgroup file, or None if it is missing or 'max'."""
    try:
        with open(path) as f:
            value = int(f.read().strip())
    except (OSError, ValueError):
        value = None
    return value

def read_cgroup_stat(path, key):
    """Return one field of a cgroup memory.stat file, or 0 if it is missing."""
    try:
        with open(path) as f:
            fields = dict(line.split() for line in f if len(line.split()) == 2)
    except OSError:
        fields = {}
    return int(fields.get(key, 0))

def read_cgroup_memory():
    """
    Return (used, limit) in bytes for the most restrictive memory cgroup that
    holds this process, or None if no limit is set. SLURM sets --mem on the job
    or step cgroup, not on the task cgroup of the process, so the cgroups above
    it are checked too. Inactive file cache can be reclaimed and is not counted
    as used.
    """
    result = None
    try:
        with open('/proc/self/cgroup') as f:
            cgroup_lines = f.read().splitlines()
    except OSError:
        cgroup_lines = []
    total_memory = psutil.virtual_memory().total
    for line in cgroup_lines:
        _, controllers, path = line.split(':', 2)
        controller = 'memory' if 'memory' in controllers.split(',') else controllers
        if result is None and controller in CGROUP_MEMORY_LAYOUTS:
            root, limit_name, usage_name, inactive_name = CGROUP_MEMORY_LAYOUTS[controller]
            directory = root + path.rstrip('/')
            while directory.startswith(root):
                limit = read_cgroup_value(os.path.join(directory, limit_name))
                usage = read_cgroup_value(os.path.join(directory, usage_name))
                # An unlimited cgroup reports 'max' (v2) or a huge number (v1)
                if (limit is not None and usage is not None and limit < total_memory
                        and (result is None or limit < result[1])):
                    inactive = read_cgroup_stat(os.path.join(directory, 'memory.stat'),
                                                inactive_name)
                    result = (usage - inactive, limit)
                directory = os.path.dirname(directory)
    return result

def memory_use_percent():
    """
    Percentage of memory in use: the larger of the use of the whole machine and,
    inside a cgroup with a memory limit (a SLURM job with --mem, for instance),
    the use of that cgroup against its limit.
    """
    percent = psutil.virtual_memory().percent
    cgroup_memory = read_cgroup_memory()
    if cgroup_memory is not None:
        used, limit = cgroup_memory
        percent = max(percent, round(100 * used / limit, 1))
    return percent

def read_shared_inputs(args):
    """
    Read the genome and parse the reference annotation once for all the chunks.
    Forked workers share them with this process instead of reading them again.
    With another start method they would be pickled for every worker, so None
    is returned and each worker reads its own, as before.
    """
    genome_dict, reference = None, None
    if get_start_method() == 'fork':
        qc_logger.info(f"Reading genome fasta {args.refFasta}")
        genome_dict = dict((r.name, r) for r in SeqIO.parse(open(args.refFasta), 'fasta'))
        # The genePred written by reference_parser is only needed while parsing
        with tempfile.TemporaryDirectory(dir=args.dir) as ref_dir:
            os.makedirs(os.path.join(ref_dir, 'logs'))
            reference = reference_parser(args.refGTF, ref_dir, args.output,
                                         list(genome_dict.keys()), args.genename,
                                         args.isoAnnotLite)
    return genome_dict, reference

def split_input_run(args, outdir, genome_dict=None, reference=None):
    SPLIT_ROOT_DIR = outdir
    if os.path.exists(SPLIT_ROOT_DIR):
        qc_logger.warning(f"{SPLIT_ROOT_DIR} directory already exists!")
    else:
        os.makedirs(SPLIT_ROOT_DIR)

    if not args.fasta:
        file_size = os.path.getsize(args.isoforms)
        target_chunk_bytes = max(1, file_size // args.chunks + (1 if file_size % args.chunks else 0))

        split_outs = []
        chunk_index = 0
        current_chunk_bytes = 0
        has_primary_feature = False
        has_non_comment_record = False
        
        d = os.path.join(SPLIT_ROOT_DIR, str(chunk_index))
        os.makedirs(d, exist_ok=True)
        current_f = open(os.path.join(d, os.path.basename(args.isoforms) + '.split' + str(chunk_index)), 'w')
        split_outs.append((os.path.abspath(d), current_f.name))

        primary_features = ('\ttranscript\t', '\tmRNA\t', '\tmatch\t', '\tcDNA_match\t')

        try:
            with open(args.isoforms, 'r') as f_in:
                for line in f_in:
                    if line.startswith('#'):
                        current_f.write(line)
                        continue
                    
                    has_non_comment_record = True
                    
                    # Check for boundary (a primary transcript line)
                    is_boundary = any(feat in line for feat in primary_features)
                    if is_boundary:
                        has_primary_feature = True
                    
                    # If we've reached the target chunk size, we only split AT a boundary
                    # so we don't sever a transcript from its downstream exons
                    if current_chunk_bytes >= target_chunk_bytes and is_boundary:
                        current_f.close()
                        chunk_index += 1
                        current_chunk_bytes = 0
                        
                        d = os.path.join(SPLIT_ROOT_DIR, str(chunk_index))
                        os.makedirs(d, exist_ok=True)
                        current_f = open(os.path.join(d, os.path.basename(args.isoforms) + '.split' + str(chunk_index)), 'w')
                        split_outs.append((os.path.abspath(d), current_f.name))

                    current_f.write(line)
                    current_chunk_bytes += len(line.encode('utf-8'))
        finally:
            if not current_f.closed:
                current_f.close()
        
        if not has_non_comment_record:
            qc_logger.error(f"Input file '{args.isoforms}' contains only comments or is empty.")
            raise ValueError(f"Input file '{args.isoforms}' contains only comments or is empty.")
            
        if not has_primary_feature:
            qc_logger.warning(f"No primary features (e.g., transcript, mRNA) detected in '{args.isoforms}'. SQANTI3 might fail to process this file correctly.")

    else:
        split_outs = split_fasta_input(args.isoforms, SPLIT_ROOT_DIR, args.chunks)
        # Align all the sequences at once, as a run without chunks does: an aligner per chunk builds its
        # own genome index, and they all run at the same time. Each chunk gets the alignments of its
        # sequences as its own SAM, so it skips the alignment. The whole SAM is removed with SPLIT_ROOT_DIR
        _, all_sam, _, _, _ = get_corr_filenames(SPLIT_ROOT_DIR, args.output)
        renamed_isoforms = align_isoforms(SPLIT_ROOT_DIR, all_sam, args.cpus, args.refFasta,
                                          args.isoforms, args.aligner_choice,
                                          args.gmap_index, args.refGTF)
        # Chunk of every cleaned up isoform ID: the chunks and the renamed file keep the input order
        chunk_sizes = []
        for _, chunk_file in split_outs:
            with open(chunk_file) as h:
                chunk_sizes.append(sum(line.startswith('>') for line in h))
        chunk_of_isoform = {}
        with open(renamed_isoforms) as h:
            isoform_ids = (line[1:].rstrip('\n') for line in h if line.startswith('>'))
            for i, size in enumerate(chunk_sizes):
                for _ in range(size):
                    chunk_of_isoform[next(isoform_ids)] = i
        # Every chunk gets the SAM header and the alignments of its isoforms
        split_sams = [open(get_corr_filenames(d, args.output)[1], 'w') for d, _ in split_outs]
        with open(all_sam) as h:
            for line in h:
                if line.startswith('@'):
                    for f in split_sams:
                        f.write(line)
                else:
                    split_sams[chunk_of_isoform[line.split('\t', 1)[0]]].write(line)
        for f in split_sams:
            f.close()

    pools = []
    
    # -t is the CPU budget of the whole run (sequence_correction() already
    # gives each worker's aligner cpus // chunks), and only the CPUs this
    # process may use count, such as a SLURM allocation, not every CPU of the
    # machine. A worker classifies on one CPU, so up to that many run at once
    if hasattr(os, 'sched_getaffinity'):
        available_cpus = len(os.sched_getaffinity(0))
    else:
        available_cpus = psutil.cpu_count(logical=True) or 1
    max_workers = max(1, min(args.cpus, available_cpus))

    pending_tasks = list(enumerate(split_outs))
    active_processes = []

    while pending_tasks or active_processes:
        # Keep only alive processes
        active_processes = [p for p in active_processes if p.is_alive()]
        
        # Check memory, against the job's limit when there is one. We need to
        # reserve some overhead (e.g., >15% free memory) for each new chunk
        mem_percent = memory_use_percent()
        can_launch_mem = mem_percent < 85  or len(active_processes) == 0
        
        # max_workers is at least 1, so a worker can always start when none is running
        can_launch_cpu = len(active_processes) < max_workers

        # Launch next chunk if resources are safe
        if pending_tasks and can_launch_cpu and can_launch_mem:
            i, (d, x) = pending_tasks.pop(0)
            qc_logger.info(f"Launching worker on {x} (Active workers: {len(active_processes)+1}, Mem usage: {mem_percent}%)")
            
            args2 = copy.deepcopy(args)
            args2.isoforms = x
            args2.novel_gene_prefix = str(i)
            args2.dir = d
            args2.report = 'skip'
            
            p = Process(target=run, args=(args2, genome_dict, reference))
            p.start()
            
            active_processes.append(p)
            pools.append(p)
            
            # Brief pause to let memory initialize before spawning another chunk.
            # Only a worker that reads its own genome and reference needs it;
            # forked workers share those of this process
            if genome_dict is None:
                time.sleep(2)
        else:
            # Wait for resources to free up: until a worker ends, checking
            # again at least every second
            wait([p.sentinel for p in active_processes], timeout=1)

    for p in pools:
        p.join()
    return [d for (d,x) in split_outs]

def concat_tables(in_files, out_file):
    """
    Concatenate tab-separated files that share a one-line header, keeping the header once.
    Missing input files are skipped. out_file is only written if at least one input exists.

    :return: number of data lines written
    """
    in_files = [f for f in in_files if os.path.exists(f)]
    if not in_files:
        return 0
    n_lines = 0
    with open(out_file, 'w') as out:
        for i, in_file in enumerate(in_files):
            with open(in_file) as h:
                header = h.readline()
                if i == 0:
                    out.write(header)
                for line in h:
                    out.write(line)
                    n_lines += 1
    return n_lines

def concat_sams(in_files, out_file):
    """
    Concatenate the SAM files of the chunks. The header (@ lines) is taken from the first chunk:
    all chunks are aligned to the same genome, so their @SQ lines are identical.
    """
    with open(out_file, 'w') as out:
        for i, in_file in enumerate(in_files):
            with open(in_file) as h:
                for line in h:
                    if line.startswith('@') and i > 0:
                        continue
                    out.write(line)

def combine_alignment_outputs(args, split_dirs):
    """
    Combine the outputs of the alignment step of each chunk (only produced with FASTA/FASTQ input):
    corrected SAM, indels table, supplementary alignments report and unmapped isoforms list.
    Without this, they are lost when the split directories are removed.
    """
    _, corrSAM, _, _, _ = get_corr_filenames(args.dir, args.output)
    split_sams = [get_corr_filenames(d, args.output)[1] for d in split_dirs]
    if not all(os.path.exists(sam) for sam in split_sams):
        qc_logger.warning("Not all chunks have an alignment SAM file. The corrected SAM will not be combined.")
        return

    concat_sams(split_sams, corrSAM)

    def indels_name(sam):
        # Same naming as calc_indels_from_sam: <corrected SAM without extension>_indels.txt
        return os.path.splitext(sam)[0] + "_indels.txt"
    concat_tables([indels_name(sam) for sam in split_sams], indels_name(corrSAM))

    # Same warnings as a run without chunks, now with the final paths (the chunks only log a debug line)
    supplementary_file = get_supplementary_name(args.dir, args.output)
    n_split = 0
    if concat_tables([get_supplementary_name(d, args.output) for d in split_dirs], supplementary_file) > 0:
        with open(supplementary_file) as h:
            n_split = len({line.split('\t')[0] for line in h} - {"isoform"})
    warn_supplementary(n_split, supplementary_file)

    unmapped_file = get_unmapped_name(args.dir, args.output)
    n_unmapped = concat_tables([get_unmapped_name(d, args.output) for d in split_dirs], unmapped_file)
    warn_unmapped(n_unmapped, unmapped_file)

def combine_split_runs(args, split_dirs, genome_dict=None):
    """
    Combine .faa, .fasta, .gtf, .classification.txt, .junctions.txt and, with FASTA/FASTQ input,
    the alignment outputs (.sam, indels, supplementary alignments, unmapped isoforms).
    Then write out the PDF report
    """
    corrGTF, _, corrFASTA, corrORF , corrCDS_GTF_GFF = get_corr_filenames(args.dir, args.output)
    outputClassPath, outputJuncPath = get_class_junc_filenames(args.dir,args.output)

    if args.isoform_hits:
        isoform_hits = open(get_isoform_hits_name(args.dir, args.output)+'_tmp', 'w')
    if args.include_ORF:
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
        if args.include_ORF:
            with open(_orf) as h: f_faa.write(h.read())
        with open(_gtf) as h: f_gtf.write(h.read())
        with open(_fasta) as h: f_fasta.write(h.read())
        with open(f"{_junc}_tmp") as h:
            if i == 0:
                f_junc_temp.write(h.readline())
            else:
                h.readline()
            f_junc_temp.write(h.read())

        if args.isoform_hits:
            _iso_hits = f"{get_isoform_hits_name(split_d, args.output)}_tmp"
            with open(_iso_hits) as h:
                if i == 0:
                    isoform_hits.write(h.readline())
                else: # skip the header line
                    h.readline()
                isoform_hits.write(h.read())
        # Retrieving the isoform object and the junctions header
        with open(_info, 'rb') as h:
            isoforms_info.update(pickle.load(h))
            headers.append(pickle.load(h))
        if args.include_ORF:
            src_dir = os.path.join(split_d, "TD2")
            dst_dir = os.path.join(args.dir, "TD2", f"TD2_{i}")
            # If the directory is not removed, further runs will create nested directories if moving directly
            if os.path.exists(dst_dir):
                shutil.rmtree(dst_dir)

            # Now move the directory.
            shutil.move(src_dir, dst_dir)
    f_fasta.close()
    f_gtf.close()
    f_junc_temp.close()

    if args.fasta:
        combine_alignment_outputs(args, split_dirs)

    # Fix novel genes and classify FSM
    isoforms_info = rename_novel_genes(isoforms_info, args.novel_gene_prefix)
    isoforms_info = classify_fsm(isoforms_info)
    
    # Process expression data natively across all processed chunks
    isoforms_info = isoform_expression_info(isoforms_info, args.expression, args.short_reads, args.dir, corrFASTA, args.cpus)
    ## FL count file
    if args.fl_count:
        isoforms_info = full_length_quantification(args.fl_count, isoforms_info, args.counts_design,
                                                   args.min_expression)
    isoforms_info,RTS_info = process_rts(isoforms_info,outputJuncPath,args.refFasta,genome_dict)

    fields_junc_cur = headers[0]
    if args.include_ORF:
        write_collapsed_GFF_with_CDS(isoforms_info, corrGTF, corrCDS_GTF_GFF)
    write_classification_output(isoforms_info, outputClassPath)
    write_junction_output(outputJuncPath, RTS_info, fields_junc_cur)
    #write omitted isoforms if requested minimum reference length is more than 0
    isoforms_info = write_omitted_isoforms(isoforms_info, args.dir, args.output,
                                            args.min_ref_len, args.is_fusion)
    #isoform hits to file if requested
    if args.isoform_hits:
        write_isoform_hits(args.dir, args.output, isoforms_info)
        isoform_hits.close()

    if args.include_ORF:
        f_faa.close()
    cleanup(outputClassPath, outputJuncPath)

    if args.report != 'skip':
        generate_report(args.saturation,args.report, outputClassPath, outputJuncPath)

# Run TUSCO benchmarking report if requested
    if hasattr(args, 'tusco') and args.tusco:
        # Pass both sample GTF (corrected) and reference GTF for IGV-like plots
        generate_tusco_report(args.tusco, outputClassPath, corrGTF, args.refGTF)
