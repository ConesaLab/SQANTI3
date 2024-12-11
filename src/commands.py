import os
import sys
import shutil

from .parsers import STARcov_parser
from .utilities.short_reads import star

utilitiesPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "utilities")
sys.path.insert(0, utilitiesPath)

GMAP_CMD = "gmap --cross-species -n 1 --max-intronlength-middle=2000000 --max-intronlength-ends=2000000 -L 3000000 -f samse -t {cpus} -D {dir} -d {name} -z {sense} {i} > {o}"
MINIMAP2_CMD = "minimap2 -ax splice --secondary=no -C5 -u{sense} -t {cpus} {g} {i} > {o}"
DESALT_CMD = "deSALT aln {dir} {i} -t {cpus} -x ccs -o {o}"
ULTRA_CMD = "uLTRA pipeline {g} {a} {i} {o_dir} --t {cpus} --prefix {prefix} --isoseq"

# To dynamically get the utilities path
def get_gmst_prog(utilities_path):
    return os.path.join(utilities_path, "gmst", "gmst.pl")


GMST_CMD = f"perl {get_gmst_prog(utilitiesPath)} -faa --strand direct --fnn --output {{o}} {{i}}"

# GTF
GTF2GENEPRED_PROG = os.path.join(utilitiesPath,"gtfToGenePred")
GFFREAD_PROG = "gffread"

if shutil.which(GTF2GENEPRED_PROG) is None:
    print("Cannot find executable {0}. Abort!".format(GTF2GENEPRED_PROG), file=sys.stderr)
    sys.exit(1)
if shutil.which(GFFREAD_PROG) is None:
    print("Cannot find executable {0}. Abort!".format(GFFREAD_PROG), file=sys.stderr)
    sys.exit(1)


# Rscript
RSCRIPTPATH = shutil.which('Rscript')
if os.system(RSCRIPTPATH + " --version")!=0:
    print("Rscript executable not found! Abort!", file=sys.stderr)
    sys.exit(1)

RSCRIPT_REPORT = '/report_qc/SQANTI3_report.R'

ISOANNOT_PROG =  os.path.join(utilitiesPath, "IsoAnnotLite_SQ3.py")


def short_reads_mapping(args):
    """
    Runs STAR for mapping short reads to the transcriptome
    """
    if args.short_reads is not None:
        print("**** Running STAR for calculating Short-Read Coverage.", file=sys.stdout)
        star_out, star_index = star(args.genome, args.short_reads, args.dir, args.cpus)
        SJcovNames, SJcovInfo = STARcov_parser(star_out)
    else:
        star_out, star_index, SJcovNames, SJcovInfo = None, None, None, None
    return(star_out, star_index, SJcovNames, SJcovInfo)
