__version__ = "0.1.0"

import sys
import os
import click
import re
import logging
from click_option_group import optgroup
from typing import Optional
import logging

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
logfile = "gms.log"
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
# gc_out = meta_out.".feature" # what?
verbose = False


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
    show_default=True,
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
    show_default=True,
)

# # Run options:
@optgroup.group("Run options")
@optgroup.option(
    "--bins",
    type=click.Choice(["0", "1", "2", "3"]),
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
    type=click.Choice(["direct", "reverse", "both"]),
    help="sequence strand to predict genes in",
    default="both",
    show_default=True,
)
@optgroup.option(
    "--order",
    type=click.IntRange(1, 100, clamp=True),
    help="markov chain order",
    default=4,
    show_default=True,
)
@optgroup.option(
    "--order_non",
    type=int,
    help="order for non-coding parameters",
    default=2,
    show_default=True,
)
@optgroup.option(
    "--gcode",
    type=click.Choice(["11", "4", "1"]),
    help="genetic code.  See https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for codes.  Currently, only tables 1, 4, and 11 are allowed.",
    default=1,
    show_default=True,
)
@optgroup.option(
    "--motif",
    type=click.Choice(["0", "1"]),
    help="iterative search for a sequence motif associated with CDS start",
    show_default=True,
    default=1,
)
@optgroup.option(
    "--width",
    type=click.IntRange(3, 100),
    help="motif width",
    show_default=True,
    default=1,
)
@optgroup.option(
    "--prestart",
    type=click.IntRange(0, 100),
    help="length of sequence upstream of translation initiation site that presumably includes the motif",
    show_default=True,
    default=6,
)
@optgroup.option(
    "--fixmotif",
    type=bool,
    is_flag=True,
    help="the motif is located at a fixed position with regard to the start motif could overlap start codon. if this option is on, it changes the meaning of the --prestart option which in this case will define the distance from start codon to motif start",
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
    help="custom parameters for GeneMarkS (default is selected based on gcode value: 'par_<gcode>.default' )",
    show_default=True,
)
@optgroup.option(
    "--gibbs",
    type=click.Choice(["1", "3"]),
    default=3,
    help="version of Gibbs sampler software (default: 3 supported versions: 1 and 3 )",
    show_default=True,
)
@optgroup.option("--test", is_flag=True, default=False, help="installation test")
@optgroup.option(
    "--identity",
    type=click.FloatRange(min=0, max=1, clamp=False),
    default=0.99,
    help="identity level assigned for termination of iterations",
    show_default=True,
)
@optgroup.option(
    "--maxitr",
    type=int,
    help="maximum number of iterations (default: 10 supported in range: >= 1)",
    show_default=True,
    default=10,
)
@optgroup.option("--verbose", is_flag=True, default=False, show_default=True)
@optgroup.option("--version", is_flag=True, default=False, show_default=True)
@click.argument("input", type=click.File("rb"))
@click.help_option(show_default=True)
@Log(verbose)
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
        order_non = 2
        offover = "0"
        gcode = "11"
        fixmotif = "0"
        prestart = 40
        width = 6
        fixmotif = 0
    if version:
        print(f"{__version__}")
    if verbose:
        # TODO: add actual logging
        # log = logging.basicConfig(filename='gmst.log', level=logging.INFO)
        pass
    pass


def cluster(feature_f, clusters):  # $gc_out, $bins
    gc_hash = dict()
    cut_off_points = []
    num_of_seq = 0
    total_length = 0
    header_to_cod_GC = dict()

    with open(feature_f, "r") as GC:
        # read in probuild output, line by line.  Should be fasta input.
        for line in GC:
            # if the line is a fasta header in the form of '>(Reference sequence name)\t(number) (number)
            text = re.match(pattern="^>(.*?)\t(\d+)\s+(\d+)", string=line)
            if text:
                header = text.group(1)  # Reference name
                length = text.group(2)  # length of sequence?
                GC = text.group(3)  # not sure, GC count?
                header_re = re.match(
                    pattern="^(.*?)\t", string=line
                )  # Dont get this one - didn't we already extract just this capture group?
                if header_re:
                    header = re.match(1)
                header_to_cod_GC[header] = GC
                num_of_seq += 1
                total_length += length
                gc_hash[GC] += length

    sorted_GC = {
        _: gc_hash[_] for _ in sorted(gc_hash)
    }  # sort the gc_hash dictionary by keys
    min_GC = sorted_GC[0]
    max_GC = sorted_GC[-1]
    print(f"min_GC={min_GC} max_GC={max_GC} total_seq_length={total_length}\n")

    previous = 0
    for key in sorted_GC:
        gc_hash[key] += previous

        if previous < total_length / 3 and gc_hash[key] >= total_length / 3:
            one_third = key
        if previous < total_length / 3 * 2 and gc_hash[key] >= total_length / 3 * 2:
            two_third = key
        if previous < total_length / 2 and gc_hash[key] >= total_length / 2:
            one_half = key
        previous = gc_hash[key]
        # TODO: uncomment when we have logging fixed
        # log.info(f"({one_third})->({gc_hash[one_third]})\n")
        # log.info(f"({one_half})->({gc_hash[one_half]})\n")
        # log.info(f"({two_third})->({gc_hash[two_third]})\n")

    if clusters == 0:
        # cluster number is not specified by user
        # automatically choose cluster number.
        if two_third - one_third > 3:
            clusters = 3
        else:
            clusters = 1
    if clusters == 3:
        if (
            (two_third - one_third) < 1
            or (max_GC - two_third) < 1
            or (one_third - min_GC) < 1
        ):
            # &Log( "Total number of sequences is not enough for training in 3 clusters!\n" )
            clusters = 1
        else:
            if gc_hash[one_third] > MIN_LENGTH:
                cut_off_points.extend((min_GC, one_third, two_third, max_GC))
            else:
                # &Log( "Total length of sequences is not enough for training in 3 clusters!\n" )
                clusters = 2

    if clusters == 2:
        if gc_hash[one_half] > MIN_LENGTH:
            cut_off_points.extend((min_GC, one_half, max_GC))
        else:
            # &Log( "Total length of sequences is not enough for training in 2 clusters!\n" )
            pass

    if clusters == 1:
        cut_off_points.extend((min_GC, max_GC))
    return (clusters, cut_off_points, header_to_cod_GC)


if __name__ == "__main__":
    if "--verbose" in sys.argv:
        verbose = True
    sys.exit(main())  # pragma: no cover
