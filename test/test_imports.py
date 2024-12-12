from utils import validate_import, validate_import_own

# Standard Library Imports
def test_import_argparse():
    validate_import("argparse", "argparse")

def test_import_sys():
    validate_import("sys", "sys")

def test_import_os():
    validate_import("os", "os")

def test_import_bisect():
    validate_import("bisect", "bisect")

def test_import_csv():
    validate_import("csv", "csv")

def test_import_glob():
    validate_import("glob", "glob")

def test_import_gzip():
    validate_import("gzip", "gzip")

def test_import_itertools():
    validate_import("itertools", "itertools")

def test_import_math():
    validate_import("math", "math")


# Third-Party Imports
def test_import_numpy():
    validate_import("numpy", "np")

def test_import_pandas():
    validate_import("pandas", "pd")

def test_import_biopython():
    validate_import("Bio.SeqIO", "SeqIO")

def test_import_bx_intervals():
    validate_import("bx.intervals", "Interval")
    validate_import("bx.intervals", "IntervalTree")


# Project-Specific Imports
import sys, os

# Add the src directory to sys.path
current_dir = os.path.dirname(os.path.abspath(__file__))
src_path = os.path.abspath(os.path.join(current_dir, ".."))
sys.path.insert(0, src_path)


def test_import_classification():
    validate_import_own("src.classification", "isoformClassification", sys.path[0])


def test_import_commands():
    validate_import_own("src.commands", "GTF2GENEPRED_PROG", sys.path[0])
    validate_import_own("src.commands", "RSCRIPTPATH", sys.path[0])
    validate_import_own("src.commands", "RSCRIPT_REPORT", sys.path[0])
    validate_import_own("src.commands", "utilitiesPath", sys.path[0])
    validate_import_own("src.commands", "GMAP_CMD", sys.path[0])
    validate_import_own("src.commands", "MINIMAP2_CMD", sys.path[0])
    validate_import_own("src.commands", "DESALT_CMD", sys.path[0])
    validate_import_own("src.commands", "ULTRA_CMD", sys.path[0])
    validate_import_own("src.commands", "GMST_CMD", sys.path[0])
    validate_import_own("src.commands", "GFFREAD_PROG", sys.path[0])


def test_import_config():
    validate_import_own("src.config", "EXP_KALLISTO_HEADERS", sys.path[0])
    validate_import_own("src.config", "EXP_RSEM_HEADERS", sys.path[0])
    validate_import_own("src.config", "FIELDS_CLASS", sys.path[0])
    validate_import_own("src.config", "FIELDS_JUNC", sys.path[0])
    validate_import_own("src.config", "seqid_fusion", sys.path[0])
    validate_import_own("src.config", "seqid_rex1", sys.path[0])
    validate_import_own("src.config", "seqid_rex2", sys.path[0])
    validate_import_own("src.config", "__version__", sys.path[0])
    validate_import_own("src.config", "__author__", sys.path[0])


def test_import_helpers():
    validate_import_own("src.helpers", "get_corr_filenames", sys.path[0])
    validate_import_own("src.helpers", "get_class_junc_filenames", sys.path[0])
    validate_import_own("src.helpers", "get_omitted_name", sys.path[0])
    validate_import_own("src.helpers", "correctionPlusORFpred", sys.path[0])
    validate_import_own("src.helpers", "write_collapsed_GFF_with_CDS", sys.path[0])
    validate_import_own("src.helpers", "write_junctionInfo", sys.path[0])
    validate_import_own("src.helpers", "get_isoform_hits_name", sys.path[0])


def test_import_multiprocessing():
    validate_import_own("multiprocessing", "Process", sys.path[0])


def test_import_parsers():
    validate_import_own("src.parsers", "get_fusion_component", sys.path[0])
    validate_import_own("src.parsers", "STARcov_parser", sys.path[0])
    validate_import_own("src.parsers", "reference_parser", sys.path[0])
    validate_import_own("src.parsers", "isoforms_parser", sys.path[0])
    validate_import_own("src.parsers", "FLcount_parser", sys.path[0])
    validate_import_own("src.parsers", "expression_parser", sys.path[0])


def test_import_qc_classes():
    validate_import_own("src.qc_classes", "genePredReader", sys.path[0])
    validate_import_own("src.qc_classes", "myQueryProteins", sys.path[0])
    validate_import_own("src.qc_classes", "myQueryTranscripts", sys.path[0])
    validate_import_own("src.qc_classes", "CAGEPeak", sys.path[0])
    validate_import_own("src.qc_classes", "PolyAPeak", sys.path[0])


def test_import_qc_pipeline():
    validate_import_own("src.qc_pipeline", "run", sys.path[0])


def test_import_utils():
    validate_import_own("src.utils", "find_polyA_motif", sys.path[0])
    validate_import_own("src.utils", "mergeDict", sys.path[0])
    validate_import_own("src.utils", "flatten", sys.path[0])
    validate_import_own("src.utils", "pstdev", sys.path[0])


def test_import_utilities():
    validate_import_own("src.utilities.cupcake.io.GFF", "collapseGFFReader", sys.path[0])
    validate_import_own("src.utilities.cupcake.io.GFF", "write_collapseGFF_format", sys.path[0])
    validate_import_own("src.utilities.cupcake.sequence.BED", "LazyBEDPointReader", sys.path[0])
    validate_import_own("src.utilities.cupcake.sequence.err_correct_w_genome", "err_correct", sys.path[0])
    validate_import_own("src.utilities.cupcake.sequence.sam_to_gff3", "convert_sam_to_gff3", sys.path[0])
    validate_import_own("src.utilities.cupcake.sequence.STAR", "STARJunctionReader", sys.path[0])
    validate_import_own("src.utilities.cupcake.tofu.compare_junctions", "compare_junctions", sys.path[0])
    validate_import_own("src.utilities.short_reads", "get_bam_header", sys.path[0])
    validate_import_own("src.utilities.short_reads", "get_ratio_TSS", sys.path[0])
    validate_import_own("src.utilities.short_reads", "get_TSS_bed", sys.path[0])
    validate_import_own("src.utilities.short_reads", "star", sys.path[0])

def test_import_statistics():
    validate_import("statistics", "mean")

def test_import_subprocess():
    validate_import("subprocess", "subprocess")

def test_import_re():
    validate_import("re", "re")

def test_import_shutil():
    validate_import("shutil", "shutil")

