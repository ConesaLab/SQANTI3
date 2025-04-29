import os, subprocess, sys
import pandas as pd
import numpy as np
import pybedtools
import re

from src.commands import run_command
from src.module_logging import qc_logger

try:
    from BCBio import GFF as BCBio_GFF
except ImportError:
    qc_logger.error("Unable to import BCBio! Please make sure bcbiogff is installed.")
    sys.exit(-1)

def star_mapping(index_dir, SR_fofn, output_dir, cpus):
    """
    Maps short reads to a reference genome using the STAR aligner.

    Parameters:
    index_dir (str): Path to the directory containing the STAR genome index.
    SR_fofn (str): Path to a file of filenames (FOFN) containing paths to short read files.
    output_dir (str): Directory where the output files will be saved.
    cpus (int): Number of CPU threads to use for the STAR aligner.

    Raises:
    FileNotFoundError: If the index directory does not exist.

    The function reads the FOFN file to get the list of short read files, checks if they are compressed,
    and then runs the STAR aligner for each sample. The output files are saved in the specified output directory.
    """

    if not os.path.exists(index_dir):
        qc_logger.error(f"Directory {index_dir} does not exist.")
        raise FileNotFoundError()

    mapping_dir = output_dir + '/STAR_mapping'
    with open(SR_fofn) as fofn:
        for line in fofn:
            files = [x.strip() for x in line.split(' ')]
            compressed = files[0].endswith('.gz')
            sample_name = os.path.splitext(os.path.basename(files[0]))[0].split('.')[0]
            sample_prefix = os.path.join(mapping_dir, sample_name)
            if not os.path.exists(f'{sample_prefix}Log.final.out'):
                qc_logger.info(f'Mapping for {sample_name}: in progress.')
                cmd = star_cmd(cpus,index_dir,files,sample_prefix,compressed) 
                logFile = f"{output_dir}/logs/STAR_mapping_{sample_name}.log"
                run_command(cmd,qc_logger,logFile,f'STAR mapping {sample_name}')
                qc_logger.info(f'Mapping for {sample_name}: done.')
            else:
                qc_logger.info(f'Mapping for {sample_name}: already completed.')


def star_cmd(cpus,index_dir,files,sample_prefix,compressed):
    base_command = [
            'STAR', '--runThreadN', str(cpus), '--genomeDir', index_dir,
            '--readFilesIn'] + files + [
            '--outFileNamePrefix', sample_prefix,
            '--alignSJoverhangMin', '8', '--alignSJDBoverhangMin', '1',
            '--outFilterType', 'BySJout', '--outSAMunmapped', 'Within',
            '--outFilterMultimapNmax', '20', '--outFilterMismatchNoverLmax', '0.04',
            '--outFilterMismatchNmax', '999', '--alignIntronMin', '20',
            '--alignIntronMax', '1000000', '--alignMatesGapMax', '1000000',
            '--sjdbScore', '1', '--genomeLoad', 'NoSharedMemory',
            '--outSAMtype', 'BAM', 'SortedByCoordinate', '--twopassMode', 'Basic'
        ]
    if compressed:
        base_command.extend(['--readFilesCommand', 'zcat'])
    return ' '.join(base_command)
        
def star(genome, SR_fofn, output_dir, cpus):
    fasta_genome = genome 
    index_dir = output_dir + '/STAR_index/'
    index_dir_tmp = index_dir + '/_STARtmp/'
    index_dir_o = index_dir + 'SAindex' 
    mapping_dir = output_dir + '/STAR_mapping/'
    if not os.path.exists(mapping_dir):
        os.makedirs(mapping_dir)
    if not os.path.exists(index_dir):
        os.makedirs(index_dir)
        if not os.path.exists(index_dir_o):
            qc_logger.info('** Running indexing.')
            cmd = ' '.join(['STAR', '--runThreadN', str(cpus), '--runMode', 'genomeGenerate', '--genomeDir', index_dir, '--genomeFastaFiles', fasta_genome, '--outTmpDir', index_dir_tmp])
            logFile = f"{output_dir}/logs/STAR_idex.log"
            run_command(cmd,qc_logger,logFile,"STAR genome indexing")
            qc_logger.info('Indexing done.')
    else:
        qc_logger.info('Index identified.')
    qc_logger.info('** Running mapping.')
    star_mapping(index_dir, SR_fofn, output_dir, cpus)
    return(mapping_dir, index_dir)

def kallisto_quantification(files,index,cpus, output_dir):
    r1 = files[0].strip()
    r2 = files[1].strip()
    if files[0][-3:] == '.gz':
        sample_name = os.path.splitext(files[0])[-2].split('/')[-1]
        sample_name = sample_name.split('.')[:-1][0]
    else:
        sample_name = os.path.splitext(files[0])[-2].split('/')[-1]
    out_prefix = output_dir + '/' + sample_name
    abundance_file = out_prefix + '/abundance.tsv'
    if not os.path.exists(abundance_file):
        if not os.path.exists(out_prefix):
            os.makedirs(out_prefix)
        qc_logger.info(f'** Running Kallisto quantification for {sample_name} sample')
        cmd = f'kallisto quant -i {index} -o {out_prefix} -b 100 -t {cpus} {r1} {r2}'
        logFile=os.path.normpath(os.path.join(output_dir,"..","logs",f"kallisto_{sample_name}.log"))
        run_command(cmd,qc_logger,logFile,"Kallisto quantification")
                           
    else:
        qc_logger.info(f"Kallisto quantification output {abundance_file} found. Using it.")
    return abundance_file


def kallisto(corrected_fasta, SR_fofn, output_dir, cpus):
    kallisto_output = output_dir + "/kallisto_output"
    kallisto_index = kallisto_output + '/kallisto_corrected_fasta.idx'
    expression_files = ''
    if not os.path.exists(kallisto_output):
        os.makedirs(kallisto_output)
    if not os.path.exists(kallisto_index):
        qc_logger.info(f'Running kallisto index {kallisto_index} using as reference {corrected_fasta}')
        cmd = f"kallisto index -i {kallisto_index} {corrected_fasta} --make-unique"
        logFile=f"{output_dir}/logs/kallisto_index.log"
        run_command(cmd,qc_logger,logFile,"Kallisto index") 
    with open(SR_fofn) as fofn:
        for line in fofn:
            files = line.split(' ')
            if len(files)==2:
                abundance = kallisto_quantification(files, kallisto_index, cpus, kallisto_output)
            else:
                qc_logger.error('SQANTI3 is only able to quantify isoforms using pair-end RNA-Seq data.\nPlease check that your fofn contains the path to both read files in a space-separated format.')
                sys.exit()
            if len(expression_files)==0:
                expression_files = abundance
            else:
                expression_files = expression_files + ',' + abundance
    return(expression_files)

def get_TSS_bed(corrected_gtf, chr_order):
    limit_info = dict(gff_type=["transcript"])
    out_directory=os.path.dirname(corrected_gtf)
    tmp_in= out_directory + '/coverage_inside_TSS.bed_tmp'
    tmp_out = out_directory + '/coverage_outside_TSS.bed_tmp'
    in_handle=open(corrected_gtf)
    with open(tmp_in,'w') as inside:
        with open(tmp_out,'w') as outside:
            for rec in BCBio_GFF.parse(in_handle, limit_info=limit_info, target_lines=1):
                chr=rec.id
                iso_id=rec.features[0].qualifiers["transcript_id"][0]
                loc=str(rec.features[0].location)
                loc=re.split(r'[\(\)\[\]\:]',loc)
                loc=list(filter(None,loc))
                strand=str(loc[2])
                if strand=='+':
                    start_in=int(loc[0])
                    end_in=int(loc[0])+100
                    start_out=int(loc[0])-101
                    end_out=int(loc[0])-1
                    if start_out<0:
                        start_out=0
                else:
                    start_in=int(loc[1])-100
                    end_in=int(loc[1])
                    start_out=int(loc[1])+1
                    end_out=int(loc[1])+101
                if end_out<=0 or start_in<=0: 
                    qc_logger.warning(f'{iso_id} will not be included in TSS ratio calculation since its TSS is located at the very beginning of the chromosome')
                else:
                    inside.write(chr + "\t" + str(start_in) + "\t" + str(end_in) + "\t" + iso_id + "\t0\t" + strand + "\n")
                    outside.write(chr + "\t" + str(start_out) + "\t" + str(end_out) + "\t" + iso_id + "\t0\t" + strand + "\n")
    in_handle.close()
    i = pybedtools.BedTool(tmp_in)
    o = pybedtools.BedTool(tmp_out)
    inside_sorted = out_directory + "/inside_TSS.bed" 
    outside_sorted = out_directory + "/outside_TSS.bed"
    i.sort(g=chr_order, output=inside_sorted)
    o.sort(g=chr_order, output=outside_sorted)
    [os.remove(i) for i in [tmp_in, tmp_out]]
    return(inside_sorted, outside_sorted)

def get_bam_header(bam):
    if not os.path.isfile(bam):
        raise FileNotFoundError(f"File {bam} not found")
    o_dir=os.path.dirname(bam)
    out=o_dir + "/chr_order.txt"
    if not os.path.isfile(out):
        cmd = rf"samtools view -H {bam} | grep '^@SQ' | sed 's/@SQ\tSN:\|LN://g'  > {out}"
        logFile=os.path.normpath(os.path.join(o_dir,"..","logs","samtools_header.log"))
        run_command(cmd,qc_logger,logFile,"samtools header")
    return(out)



def get_ratio_TSS(inside_bed, outside_bed, replicates, chr_order, metric): 
## calculate the average coverage per sample for in and out beds. Calculate each ratio
## get ratios across replicates and return it as a dictionary
    def process_coverage(cov, column_name):
        data = [(entry.name, float(entry[6])) for entry in cov]
        df = pd.DataFrame(data, columns=['id', column_name])
        if column_name == 'inside':
            df.loc[df[column_name] < 3, column_name] = np.nan
        return df

    qc_logger.info('BAM files identified: '+str(replicates))
    out_TSS_file = os.path.dirname(inside_bed) + "/ratio_TSS.csv"
    in_bed = pybedtools.BedTool(inside_bed)
    out_bed = pybedtools.BedTool(outside_bed)
    ratio_rep_df = None
    for b,bam_file in enumerate(replicates):
        in_cov = in_bed.coverage(bam_file, sorted=True, g=chr_order)
        out_cov = out_bed.coverage(bam_file, sorted=True, g=chr_order)
        inside_df = process_coverage(in_cov, 'inside')
        outside_df = process_coverage(out_cov, 'outside')
        merged = pd.merge(inside_df, outside_df, on="id")
        merged['ratio_TSS'] = (merged['inside'] + 0.01) / (merged['outside'] + 0.01)
        if ratio_rep_df is None:
            ratio_rep_df = merged[['id', 'ratio_TSS']]
        else:
            ratio_rep_df = pd.merge(ratio_rep_df, merged[['id', 'ratio_TSS']], on='id')
        
        ratio_rep_df = ratio_rep_df.rename(columns={'ratio_TSS': f'ratio_TSS_{b}'})
    ratio_rep_df.iloc[:, 1:] = ratio_rep_df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')
    ratio_rep_df.to_csv(out_TSS_file, index=False)

    # Convert all columns but id to numeric, coercing any non-numeric values to NaN
    if metric == "mean":
        ratio_rep_df['return_ratio'] = ratio_rep_df.mean(axis=1, numeric_only=True, skipna=True)
    elif metric == "3quartile":
        ratio_rep_df['return_ratio'] = ratio_rep_df.quantile(q=0.75, axis=1, numeric_only=True, skipna=True)
    elif metric == "max":
        ratio_rep_df['return_ratio'] = ratio_rep_df.max(axis=1, numeric_only=True, skipna=True)
    elif metric == "median":
        ratio_rep_df['return_ratio'] = ratio_rep_df.median(axis=1, numeric_only=True, skipna=True)
    else:
        raise ValueError("Invalid value for 'metric'. Use 'mean', '3quartile', 'max', or 'median'.")

    ratio_rep_df = ratio_rep_df[['id','return_ratio']]
    ratio_rep_dict = ratio_rep_df.set_index('id').T.to_dict()
    # [os.remove(i) for i in [inside_bed, outside_bed]]
    # print('Temp files removed.\n')
    return(ratio_rep_dict) 
