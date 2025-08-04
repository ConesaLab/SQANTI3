import os

from Bio import SeqIO

from src.utilities.cupcake.io.BioReaders import GMAPSAMReader
from src.utilities.cupcake.io.GFF import collapseGFFReader, write_collapseGFF_format

from src.commands import run_command
from src.module_logging import filter_logger


def filter_isoforms(filename,prefix,ids_to_keep):
    fafq_type = 'fasta'
    with open(filename) as h:
            if h.readline().startswith('@'): fafq_type = 'fastq'
    fout=open(prefix + '.filtered.' + fafq_type, 'w')
    for r in SeqIO.parse(open(filename), fafq_type):
        if r.id in ids_to_keep:
            SeqIO.write(r, fout, fafq_type)
    fout.close()
    filter_logger.info(f"Output written to: {fout.name}")

def filter_gtf(filename,prefix,ids_to_keep):
    outputGTF = prefix + '.filtered.gtf'
    with open(outputGTF, 'w') as f:
        for r in collapseGFFReader(filename):
            # Fix negative strand coordinates
            if (r.strand == '-') and (len(r.ref_exons) > 0): 
                r.ref_exons.reverse()
                r.start , r.end = r.ref_exons[0].start, r.ref_exons[-1].end # Positions get shifted by 1
            if r.seqid in ids_to_keep:
                write_collapseGFF_format(f, r)
        filter_logger.info(f"Output written to: {f.name}")
    f.close()

def filter_sam(filename,prefix,ids_to_keep):
    outputSAM = prefix + '.filtered.sam'
    with open(outputSAM, 'w') as f:
        for r in GMAPSAMReader(filename,has_header=True):
            if r.qID in ids_to_keep:
                f.write(r)
        filter_logger.info(f"Output written to: {f.name}")
    f.close()

def filter_gff3(filename,prefix,inclusion_f):
    outputGFF3 = prefix + '.filtered.gff3'
    awk_cmd = """awk 'FNR==NR {{ a[$1]; next }} ($1 in a)' {l} {g} > {o}""".format(l=inclusion_f, g=filename, o=outputGFF3)
    logFile = os.path.join(os.path.dirname(prefix), 'logs', 'filter_gff3.log')
    run_command(awk_cmd,filter_logger,logFile,"Filtering GFF3 file")
    filter_logger.info(f"Output written to: {outputGFF3}")

def filter_faa(filename,prefix,ids_to_keep):
    outputFAA = prefix + '.filtered.faa'
    with open(outputFAA, 'w') as f:
        for r in SeqIO.parse(open(filename), 'fasta'):
            if r.id in ids_to_keep:
                f.write(">{0}\n{1}\n".format(r.description, r.seq))
        filter_logger.info(f"Output written to: {f.name}")
    f.close()