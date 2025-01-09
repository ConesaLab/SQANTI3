import os,sys,copy,csv
import pandas as pd

from multiprocessing import Process
from Bio import SeqIO
from .utilities.cupcake.io.GFF import collapseGFFReader, write_collapseGFF_format

from .qc_pipeline import run
from .helpers import get_corr_filenames, get_class_junc_filenames
from .qc_output import generate_report
# TODO: Do the split based on isoform ID groups, not on pure numbers


def get_split_dir(outdir,prefix):
    split_prefix=os.path.join(os.path.abspath(outdir), prefix)
    split_directory = split_prefix+'_splits/'
    return split_directory

def split_input_run(args,outdir):

    SPLIT_ROOT_DIR = outdir
    if os.path.exists(SPLIT_ROOT_DIR):
        print("WARNING: {0} directory already exists!".format(SPLIT_ROOT_DIR), file=sys.stderr)
    else:
        os.makedirs(SPLIT_ROOT_DIR)

    if not args.fasta:
        # check if the code work if not use np to read the isoform file
        try:
            recs = [r for r in collapseGFFReader(args.isoforms)]
        except Exception as e:
            # read the file args.isoforms as a np file since the GTF file is not in a formal format
            recs_df = pd.read_csv(args.isoforms, sep='\t', comment='#', header=None)
            # Extract transcript IDs and assign them to the new column
            for i, value in enumerate(recs_df.iloc[:, 8]):
                parts = value.split('; ')
                for part in parts:
                    if 'transcript_id' in part:
                        transcript_id = part.split('"')[1]
                        recs_df.at[i, 'transcript_id'] = transcript_id
                        break  # Assuming only one transcript_id per row, break
            recs = {transcript_id: group.iloc[:, :-1] for transcript_id, group in recs_df.groupby('transcript_id')}

        n = len(recs)
        # if resc is empty, then the file is not in the correct format and ask user to check
        if n == 0:
            print("The input file is not in the correct format, please check the file contains transcript_id in "
                  "column 9 and try again")
            sys.exit(1)
        chunk_size = n // args.chunks + (n % args.chunks > 0)
        split_outs = []

        for i in range(args.chunks):
            if i * chunk_size >= n:
                break
            d = os.path.join(SPLIT_ROOT_DIR, str(i))
            try:
                os.makedirs(d)
            except FileExistsError:
                pass
            f = open(os.path.join(d, os.path.basename(args.isoforms) + '.split' + str(i)), 'w')
            if type(recs) == dict:
                for key in sorted(recs.keys())[i * chunk_size: min((i + 1) * chunk_size, n)]:
                    # Append each DataFrame to the file without index and with tab separation
                    recs[key].to_csv(f, sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE, escapechar='\\')
            else:
                for j in range(i * chunk_size, min((i + 1) * chunk_size, n)):
                        write_collapseGFF_format(f, recs[j])
            f.close()
            split_outs.append((os.path.abspath(d), f.name))
    else:
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
        print("launching worker on {0}....".format(x))
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
    f_class = open(outputClassPath, 'w')
    f_junc = open(outputJuncPath, 'w')
    f_cds_gtf_gff = open(corrCDS_GTF_GFF, 'w')

    for i,split_d in enumerate(split_dirs):
        _gtf, _, _fasta, _orf , _CDS_GTF_GFF = get_corr_filenames(split_d,args.output)
        _class, _junc = get_class_junc_filenames(split_d,args.output)
        if not args.skipORF:
            with open(_orf) as h: f_faa.write(h.read())
        with open(_gtf) as h: f_gtf.write(h.read())
        with open(_fasta) as h: f_fasta.write(h.read())
        with open(_class) as h:
            if i == 0:
                f_class.write(h.readline())
            else:
                h.readline()
            f_class.write(h.read())
        with open(_junc) as h:
            if i == 0:
                f_junc.write(h.readline())
            else:
                h.readline()
            f_junc.write(h.read())
        with open(_CDS_GTF_GFF) as h:
            if i == 0: # This if condition checks if its the first file to write the header or not in the final file
                f_cds_gtf_gff.write(h.readline())
            else:
                h.readline()
            f_cds_gtf_gff.write(h.read())

    f_fasta.close()
    f_gtf.close()
    f_class.close()
    f_junc.close()
    f_cds_gtf_gff.close()
    if not args.skipORF:
        f_faa.close()

    if args.report != 'skip':
        generate_report(args.saturation,args.report, outputClassPath, outputJuncPath)
