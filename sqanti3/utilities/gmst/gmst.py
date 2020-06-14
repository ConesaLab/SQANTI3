__version__ = "0.1.0"

import sys
import os
import click
import regex
from click_option_group import optgroup
from typing import Optional

BINS = "1|2|3|0"
SHAPE_TYPE = "linear|circular|partial"
GENETIC_CODE = "11|4|1"
STRAND = "direct|reverse|both"
MIN_HEURISTIC_GC = 30
MAX_HEURISTIC_GC = 70
OUTPUT_FORMAT = "LST|GFF"
MIN_LENGTH = 10000

# TODO: redefine these as needed within the code
minGC = 30
maxGC = 70
metaout = "meta.lst"
logfile  = "gms.log"
seq = "sequence"
start_prefix = "startseq."
gibbs_prefix = "gibbs_out."
mod_prefix = "itr_"
mod_suffix = ".mod"
hmmout_prefix = "itr_"
hmmout_suffix = ".lst"
out_name = "GeneMark_hmm.mod"
out_name_heuristic = "GeneMark_hmm_heuristic.mod"
out_suffix = "_hmm.mod"
out_suffix_heu = "_hmm_heuristic.mod"
fnn_out = ""
faa_out = ""

meta_out = "initial.meta.lst"
gc_out = $meta_out.".feature"

@click.command()
@optgroup.group("Output options", help="output is in current working directory")
@optgroup.option(
    "--output",
    "-o",
    type=str,
    help="output file with predicted gene coordinates by GeneMarh.hmm and species parameters derived by GeneMarkS-T. If not provided, the file name will be taken from the input file. GeneMark.hmm can be executed independently after finishing GeneMarkS training.This method may be the preferable option in some situations, as it provides accesses to GeneMarh.hmm options.",
    required=False,
)
@optgroup.option(
    "--format",
    "outputformat",
    type=str,
    help="output coordinates of predicted genes in LST or GTF format.",
    default="LST",
    show_default=True
)
@optgroup.option(
    "--fnn",
    is_flag=True,
    help="create file with nucleotide sequence of predicted genes",
    default=False,
    show_default=True,
)
@optgroup.option(
    "--faa",
    is_flag=True,
    help="create file with protein sequence of predicted genes",
    default=False,
    show_default=True,
)
@optgroup.option(
    "--clean",
    is_flag=True,
    help="delete all temporary files",
    default=True,
    show_default=True
)

# # Run options:
@optgroup.group("Run options")
@optgroup.option(
    "--bins",
    type=click.Choice(["0","1","2","3"]),
    help="number of clusters for inhomogeneous genome. Use 0 for automatic clustering",
    default=0,
    show_default=True,
)
@optgroup.option(
    "--filter",
    "filterseq",
    type=int,
    help="keep at most one prediction per sequence",
    show_default=True,
    default=True,
)
@optgroup.option(
    "--strand",
    type=click.Choice(["direct","reverse","both"]),
    help="sequence strand to predict genes in",
    default="both",
    show_default=True,
)
@optgroup.option(
    "--order",
    type=click.IntRange(1,100, clamp=True),
    help="markov chain order",
    default=4,
    show_default=True,
)
@optgroup.option(
    "--order_non",
    type=int,
    help="order for non-coding parameters",
    default = 2,
    show_default=True,
)
@optgroup.option(
    "--gcode",
    type=click.Choice(["11","4","1"]),
    help="genetic code.  See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for codes.  Currently, only tables 1, 4, and 11 are allowed.",
    default=1,
    show_default=True,
)
@optgroup.option(
    "--motif",
    type=click.Choice(['0','1']),
    help="iterative search for a sequence motif associated with CDS start",
    show_default=True,
    default=1
)
@optgroup.option(
    "--width", 
    type=click.IntRange(3, 100),
    help="motif width",show_default=True,default=1,
)
@optgroup.option(
    "--prestart",
    type=click.IntRange(0,100),
    help="length of sequence upstream of translation initiation site that presumably includes the motif",show_default=True,default=6
)
@optgroup.option(
    "--fixmotif",
    type=bool,
    is_flag=True,
    help="the motif is located at a fixed position with regard to the start; motif could overlap start codon. if this option is on, it changes the meaning of the --prestart option which in this case will define the distance from start codon to motif start",
    show_default=True,
    default=True,
)
@optgroup.option(
    "--offover",
    type=bool,
    is_flag=True,
    default=True,
    help="prohibits gene overlap",
    show_default=True,
)
@optgroup.group("Combined output and run options")
@optgroup.option(
    "--prok",
    is_flag=True,
    default=False,
    help="to run program on prokaryotic transcripts (this option is the same as: --bins 1 --filter 0 --order 2 --order_non 2 --gcode 11 -width 6 --prestart 40 --fixmotif 0)",
    show_default=True,
)
@optgroup.group("Test/developer options")
@optgroup.option(
    "--par",
    type=str,
    help="custom parameters for GeneMarkS (default is selected based on gcode value: 'par_<gcode>.default' )",show_default=True,
)
@optgroup.option(
    "--gibbs",
    type=click.Choice(['1','3']),
    default=3,
    help="version of Gibbs sampler software (default: $gibbs_version; supported versions: 1 and 3 )" ,show_default=True,
)
@optgroup.option("--test", is_flag=True, default=False, help="installation test")
@optgroup.option(
    "--identity",
    type=click.FloatRange(min = 0, max = 1, clamp = False),
    default=0.99,
    help="identity level assigned for termination of iterations",show_default=True,
)
@optgroup.option(
    "--maxitr",
    type=int,
    help="maximum number of iterations (default: $maxitr; supported in range: >= 1)",show_default=True,default=10,
)
@optgroup.option("--verbose", is_flag=True, default=False, show_default=True,)
@optgroup.option("--version", is_flag=True, default=False,show_default=True,)
@click.argument('input', type=click.File('rb'))
@click.help_option(show_default=True)
def main(
    sequence_file_name: str,
    outputformat: Optional[str] = None,
    format: str = "LST",
    fnn: bool = False,
    faa: bool = False,
    clean: bool = True,
    bins: int = 0,
    prok: bool = False,
    filterseq: int = 1,
    strand: str = "both",
    order: int = 4,
    order_non: int = 2,
    gcode: int = 1,
    motif: int = 1,
    width: int = 12,
    prestart: int = 6,
    fixmotif: int = 1,
    offover: int = 1,
    par: Optional[str] = None,
    gibbs: int = 3,
    test: bool = False,
    identity: float = 0.99,
    maxitr: int = 10,
    verbose: bool = False,
    version: bool = False,
) -> None:
    bins = int(bins)

    if output is None:
        base = os.path.basename(input)
        output = os.path.splitext(base)[0]
    if par:
        bins = 1
        filter = 0
        order = 2
        order_non=2
        offover = '0'
        gcode="11"
        fixmotif="0"
        prestart=40
        width=6
        fixmotif=0
    if version:
        print(f"{__version__}")
    
    pass


def cluster(feature_f, clusters):  # $gc_out, $bins
    gc_hash = dict()
    cut_off_points = []
    # not sure yet my ($min_GC, $max_GC, $one_third, $two_third, $one_half); # is this a tuple?
    num_of_seq = 0
    total_length = 0
    header_to_cod_GC = dict()

    with open(feature_f, 'r') as GC:
        while (<GC>):
            if ($_ !~ /^>(.*?)\t(\d+)\s+(\d+)/);
                my $header = $1;
                my $length = $2;
                my $GC = $3;
            elif($header =~ /^(.*?)\t/):
                $header = $1;
            $header_to_cod_GC{$header} = $GC;
            $num_of_seq ++;
            $total_length += $length;
            $gc_hash{$GC} += $length;

def check_args(func):
    def wrapper(*args, **kwargs):


def accepts(*types):
    def check_accepts(f):
        assert outputformat in ("LST", "GTF"), print(f"Error: format {outputformat} is not supported")
        # if type(order) not int and len(order) > 2:
        #     print(f"Error: Markov chain order {order} is not supported")
        #     exit
        # if gcode not in ("1", "4", "11"):
        #     print(f"Error: genetic code {gcode} is not supported")
        #     exit
        # if motif not in (0, 1):
        #     print( "Error: in value of motif option")
        #     exit
        # if clean not in (0, 1):
        #     print("Error: in value of clean option")
        #     exit
        # if filterseq not in (0, 1):
        #     print("Error: in value of filter option")
        #     exit
        # if fixmotif not in (0, 1):
        #     print("Error: in value of fixmotif option")
        #     exit
        # if offover not in (0, 1):
        #     print("Error: in value of offover option")
        #     exit
        # if bins not in ("0","1","2","3"):
        #     print("Error: in value of bin option")
        #     exit
        # if width !~ /^\d+$/ ):
        #     print("Error: in value of motif width option")
        #     exit
        # if prestart !~ /^\d+$/ ):
        #     print("Error: in value of prestart option")
        #     exit
        # if (identity !~ /^\d?\.?\d+$/ )||identity > 1 )):
        #     print("Error: in value of identity option")
        #     exit
        # if maxitr !~ /^\d+$/ ):
        #     print("Error: in value of maximum number of itteration")
        #     exit
        # if STRAND !~ /\b$strand\b/ ):
        #     print("Error: strand [$strand] is not supported")
        #     exit
        # if gibbs_version != 1)&&($gibbs_version != 3))
        #     print("Error: in specification of Gibbs version")
        #     exit
        new_f.__name__ = f.__name__
        return new_f
    return check_accepts

if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
