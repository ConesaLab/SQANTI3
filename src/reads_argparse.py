import argparse
from src.config import __reads_version__

def reads_argparser():
    #arguments
    parser = argparse.ArgumentParser(description="Structural and Quality Annotation of Novel Transcript Isoforms")

    pr = parser.add_argument_group("Required arguments")
    pr.add_argument('--refFasta',  required=True, help='Reference genome (Fasta format)')
    pr.add_argument('--refGTF',  required= True, help='Reference annotation file (GTF format)')
    pr.add_argument('-de', '--design', type=str, dest="inDESIGN" ,required=True, help='Path to design file, must have sampleID and file_acc column.')
    
    pio = parser.add_argument_group("Input/Output options")
    pio.add_argument('-i', '--raw_data_dir', type=str, dest='input_dir', default = './', help = '\t\tPath to directory where fastq/GTF files are stored (for running SQANTI3 from scratch). Default: Directory where the script was run.')
    pio.add_argument('-p','--prefix', type=str, dest="PREFIX", required=False, help='SQANTI-reads output filename prefix. Default: sqantiReads')
    pio.add_argument('-d', '--sqanti_dirs', type=str, dest='sqanti_dirs', help='\t\tDirectory containing existing SQANTI3 output folders. Use this to skip re-running QC and proceed directly to aggregation.', default = None, required=False)
    pio.add_argument('-o','--output', type=str, dest="OUTPUT", help='\t\tDirectory for output sqanti_reads files (plots, tables, design file). Default: Directory where the script was run.', default = "./", required=False)
    pio.add_argument('--report', type=str, choices = ["pdf", "html", "both"], default = 'pdf', help = "\t\tDefault: pdf")
    pio.add_argument('--all_tables', dest="ALLTABLES", action='store_true', help='Export all output tables. Default tables are gene counts, ujc counts, length_summary, cv and underannotated gene tables')
    pio.add_argument('--pca_tables', dest="PCATABLES", action='store_true', help='Export table for making PCA plots')

    pf = parser.add_argument_group("Filtering options")
    pf.add_argument('--min_ref_len', type=int, default=0, help="\t\tMinimum reference transcript length. Default: 0 bp")
    pf.add_argument('-ge','--gene_expression', type=int, dest="ANNOTEXP", required=False, help='Expression cut off level for determining underannotated genes. Default = 100', default = 100)
    pf.add_argument('-je','--jxn_expression', type=int, dest="JXNEXP", required=False, help='Coverage threshold for detected reference donors and acceptor. Default = 10', default = 10)
    pf.add_argument('-pc','--perc_coverage', type=int, dest="PERCCOV", required=False, help='Percent gene coverage of UJC for determining well-covered unannotated transcripts. Default = 20', default = 20)
    pf.add_argument('-pj','--perc_junctions', type=int, dest="PERCMAXJXN", required=False, help='Percent of the max junctions in gene for determining near full-length putative novel transcripts. Default = 80', default = 80)
   
    pa = parser.add_argument_group("Analysis options")
    pa.add_argument('--aligner_choice', type=str, choices=['minimap2', "uLTRA"], default='minimap2', help="\t\tDefault: minimap2")
    pa.add_argument('-s','--sites', type=str, default="ATAC,GCAG,GTAG", help='\t\tSet of splice sites to be considered as canonical (comma-separated list of splice sites). Default: GTAG,GCAG,ATAC.', required=False)
    pa.add_argument('--skip_hash', dest="SKIPHASH", action='store_true', help='Skip the hashing step')

    pv = parser.add_argument_group("Visualization options")
    pv.add_argument('-f', '--factor', type=str, dest="inFACTOR" ,required=False, help='This is the column name that plots are to be faceted by. Default: None')
    pv.add_argument('-fl','--factor_level', type=str, dest="FACTORLVL", required=False, help='Factor level to evaluate for underannotation', default = None)

    pp = parser.add_argument_group("Performance options")
    pp.add_argument('-t', '--cpus', default=10, type=int, help='\t\tNumber of threads used during alignment by aligners. Default: 10')
    pp.add_argument('-n', '--chunks', default=1, type=int, help='\t\tNumber of chunks to split SQANTI3 analysis in for speed up. Default: 1')

    # Check because this might not be needed
    parser.add_argument('--force_id_ignore', action="store_true", default=False, help="\t\t Allow the usage of transcript IDs non related with PacBio's nomenclature (PB.X.Y)")


    poa = parser.add_argument_group("Optional arguments")
    poa.add_argument('--verbose', help = 'If verbose is run, it will print all steps, by default it is FALSE', action="store_true")
    poa.add_argument('-v', '--version', help="Display program version number.", action='version', version='sqanti-reads '+str(__reads_version__))


    return parser