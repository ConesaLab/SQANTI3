import os, subprocess, sys
import pandas
import numpy as np
import pybedtools
import re
try:
    from BCBio import GFF as BCBio_GFF
except ImportError:
    print("Unable to import BCBio! Please make sure bcbiogff is installed.", file=sys.stderr)
    sys.exit(-1)

def star_mapping(index_dir, SR_fofn, output_dir, cpus): #added cpus argument for num of threads
    mapping_dir = output_dir + '/STAR_mapping'
    with open(SR_fofn) as fofn:
        for line in fofn:
            files = [x.strip() for x in line.split(' ')]
            compressed = False
            if files[0][-3:] == '.gz':
                compressed = True
            if compressed :
                sample_name = os.path.splitext(files[0])[-2].split('/')[-1]
                sample_name = sample_name.split('.')[:-1][0]
            else:
                sample_name = os.path.splitext(files[0])[-2].split('/')[-1]
            sample_prefix = mapping_dir + '/' + sample_name
            if not os.path.exists(sample_prefix + 'Log.final.out'):
                print('Mapping for ', sample_name, ': in progress...')
                if not compressed:
                    if len(files) == 1:
                        print('Mapping for ', sample_name, ': done.')
                        subprocess.call(['STAR', '--runThreadN',  str(cpus), '--genomeDir', index_dir, '--readFilesIn', files[0], '--outFileNamePrefix', sample_prefix,'--alignSJoverhangMin', '8', '--alignSJDBoverhangMin', '1', '--outFilterType', 'BySJout', '--outSAMunmapped', 'Within', '--outFilterMultimapNmax', '20', '--outFilterMismatchNoverLmax', '0.04', '--outFilterMismatchNmax', '999', '--alignIntronMin', '20', '--alignIntronMax', '1000000', '--alignMatesGapMax', '1000000', '--sjdbScore', '1', '--genomeLoad', 'NoSharedMemory', '--outSAMtype', 'BAM', 'SortedByCoordinate', '--twopassMode', 'Basic'])
                    else:
                        print('Mapping for ', sample_name, ': done.')
                        subprocess.call(['STAR', '--runThreadN', str(cpus),  '--genomeDir', index_dir, '--readFilesIn', files[0], files[1], '--outFileNamePrefix', sample_prefix,'--alignSJoverhangMin', '8', '--alignSJDBoverhangMin', '1', '--outFilterType', 'BySJout', '--outSAMunmapped', 'Within', '--outFilterMultimapNmax', '20', '--outFilterMismatchNoverLmax', '0.04', '--outFilterMismatchNmax', '999', '--alignIntronMin', '20', '--alignIntronMax', '1000000', '--alignMatesGapMax', '1000000', '--sjdbScore', '1', '--genomeLoad', 'NoSharedMemory', '--outSAMtype', 'BAM', 'SortedByCoordinate', '--twopassMode', 'Basic'])
                else:
                    if len(files) == 1:
                        print('Mapping for ', sample_name, ': done.')
                        subprocess.call(['STAR', '--runThreadN',  str(cpus), '--genomeDir', index_dir, '--readFilesIn', files[0], '--outFileNamePrefix', sample_prefix,'--alignSJoverhangMin', '8', '--alignSJDBoverhangMin', '1', '--outFilterType', 'BySJout', '--outSAMunmapped', 'Within', '--outFilterMultimapNmax', '20', '--outFilterMismatchNoverLmax', '0.04', '--outFilterMismatchNmax', '999', '--alignIntronMin', '20', '--alignIntronMax', '1000000', '--alignMatesGapMax', '1000000', '--sjdbScore', '1', '--genomeLoad', 'NoSharedMemory', '--outSAMtype', 'BAM', 'SortedByCoordinate', '--readFilesCommand', 'zcat', '--twopassMode', 'Basic'])
                    else:
                        print('Mapping for ', sample_name, ': done.')
                        subprocess.call(['STAR', '--runThreadN', str(cpus), '--genomeDir', index_dir, '--readFilesIn', files[0], files[1], '--outFileNamePrefix', sample_prefix,'--alignSJoverhangMin', '8', '--alignSJDBoverhangMin', '1', '--outFilterType', 'BySJout', '--outSAMunmapped', 'Within', '--outFilterMultimapNmax', '20', '--outFilterMismatchNoverLmax', '0.04', '--outFilterMismatchNmax', '999', '--alignIntronMin', '20', '--alignIntronMax', '1000000', '--alignMatesGapMax', '1000000', '--sjdbScore', '1', '--genomeLoad', 'NoSharedMemory', '--outSAMtype', 'BAM', 'SortedByCoordinate', '--readFilesCommand', 'zcat', '--twopassMode', 'Basic'])


def star(genome, SR_fofn, output_dir, cpus):
    fasta_genome = genome #Fasta Format already checked
    index_dir = output_dir + '/STAR_index/'
    index_dir_tmp = index_dir + '/_STARtmp/'
    index_dir_o = index_dir + 'SAindex' 
    mapping_dir = output_dir + '/STAR_mapping/'
    print('START running STAR...')
    if not os.path.exists(mapping_dir):
        os.makedirs(mapping_dir)
    if not os.path.exists(index_dir):
        os.makedirs(index_dir)
        if not os.path.exists(index_dir_o):
            print('Running indexing...')
            subprocess.call(['STAR', '--runThreadN', str(cpus), '--runMode', 'genomeGenerate', '--genomeDir', index_dir, '--genomeFastaFiles', fasta_genome, '--outTmpDir', index_dir_tmp])
            print('Indexing done.')
    else:
        print('Index identified. Proceeding to mapping.')
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
        print('Running Kallisto quantification for {0} sample'.format(sample_name))
        os.system('kallisto quant -i ' + index + ' -o ' + out_prefix + ' -b 100 -t ' + str(cpus) + ' ' + r1 + ' ' + r2)
        #subprocess.call(['kallisto', 'quant', '-i', index, '-o', out_prefix, '-b', '100', '-t', str(cpus), r1, r2 ],shell=True)
    else:
        print("Kallisto quantification output {0} found. Using it...".format(abundance_file))
    return(abundance_file)


def kallisto(corrected_fasta, SR_fofn, output_dir, cpus):
    kallisto_output = output_dir + "/kallisto_output"
    kallisto_index = kallisto_output + '/kallisto_corrected_fasta.idx'
    expression_files = ''
    if not os.path.exists(kallisto_output):
        os.makedirs(kallisto_output)
    if not os.path.exists(kallisto_index):
        print('Running kallisto index {0} using as reference {1}'.format(kallisto_index, corrected_fasta))
        os.system('kallisto index -i ' + kallisto_index + ' ' + corrected_fasta + ' --make-unique')
        #subprocess.call(['kallisto', 'index', '-i', kallisto_index, corrected_fasta, '--make-unique'],shell=True)
    with open(SR_fofn) as fofn:
        for line in fofn:
            files = line.split(' ')
            if len(files)==2:
                abundance = kallisto_quantification(files, kallisto_index, cpus, kallisto_output)
            else:
                print('SQANTI3 is only able to quantify isoforms using pair-end RNA-Seq data. Please check that your fofn contains the path to both read files in a space-separated format.')
                sys.exit()
            if len(expression_files)==0:
                expression_files = abundance
            else:
                expression_files = expression_files + ',' + abundance
    return(expression_files)


#if args.SR_fofn and not args.coverage:
#	import utilities/short_reads as sr
#	print ('Short-read files provided.')
#	args.coverage = sr.star(args.genome, args.SR_fofn, args.dir)
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
                loc=re.split('[\(\)\[\]\:]',loc)
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
                    print('{iso} will not be included in TSS ratio calculation since its TSS is located at the very beginning of the chromosome'.format(iso=iso_id))
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
    os.system("rm {i} {o}".format(i=tmp_in , o=tmp_out))
    return(inside_sorted, outside_sorted)

def get_bam_header(bam):
    o_dir=os.path.dirname(bam)
    out=o_dir + "/chr_order.txt"
    os.system("samtools view -H {b} | grep '^@SQ' | sed 's/@SQ\tSN:\|LN://g'  > {o}".format(b=bam, o=out))
    return(out)

def get_ratio_TSS(inside_bed, outside_bed, replicates, chr_order, metric): 
## the idea would be to first calculate the average coverage per sample for in and out beds. Calculate each ratio
## get the maximum the ratios across replicates and return it as a dictionary
    print('BAM files identified: '+str(replicates))
    out_TSS_file = os.path.dirname(inside_bed) + "/ratio_TSS.csv"
    in_bed = pybedtools.BedTool(inside_bed)
    out_bed = pybedtools.BedTool(outside_bed)
    for b in [*range(0,len(replicates))]:
        bam_file=replicates[b] 
        in_cov = in_bed.coverage(bam_file, sorted=True, g=chr_order)
        out_cov = out_bed.coverage(bam_file, sorted=True, g=chr_order)
        inside_df = pandas.DataFrame(columns=['id','inside'])
        for entry in in_cov:
            new_entry = pandas.DataFrame({'id' : [entry.name] , 'inside' : [float(entry[6])]})
            if (new_entry['inside'] < 3).bool():
                new_entry['inside'] = np.nan
            inside_df = pandas.concat([inside_df,new_entry], ignore_index=True)
        outside_df = pandas.DataFrame(columns=['id','outside'])
        for entry in out_cov:
            new_entry = pandas.DataFrame({'id' : [entry.name] , 'outside' : [float(entry[6])]})
            outside_df = pandas.concat([outside_df, new_entry], ignore_index=True)
        merged = pandas.merge(inside_df, outside_df, on="id")
        merged['ratio_TSS'] = (merged['inside']+0.01)/(merged['outside']+0.01)
        merged['ratio_TSS'] = pandas.to_numeric(merged['ratio_TSS'])
        if b == 0 :
            ratio_rep_df = merged['id']
        ratio_rep_df = pandas.merge(ratio_rep_df, merged[['id','ratio_TSS']], on='id')
        renamed_ratioTSS = "ratio_TSS_" + str(b)
        ratio_rep_df = ratio_rep_df.rename(columns={'ratio_TSS':renamed_ratioTSS})

    #ratio_rep_df.to_csv(out_TSS_file, index=False)
    
    # use metric value to get the final ratio_TSS that will be recorded in the classification file
    if metric == "mean":
        ratio_rep_df['return_ratio'] = ratio_rep_df.mean(axis=1, numeric_only=True, skipna=True)
    elif metric == "3quartile":
        dratio_rep_df['return_ratio'] = ratio_rep_df.quantile(q=0.75, axis=1, numeric_only=True, skipna=True)
    elif metric == "max":
        ratio_rep_df['return_ratio'] = ratio_rep_df.max(axis=1, numeric_only=True, skipna=True)
    elif metric == "median":
        ratio_rep_df['return_ratio'] = ratio_rep_df.median(axis=1, numeric_only=True, skipna=True)
    else:
        raise ValueError("Invalid value for 'metric'. Use 'mean', '3quartile', 'max', or 'median'.")

    #ratio_rep_df['max_ratio_TSS']=ratio_rep_df.max(axis=1, numeric_only=True)
    ratio_rep_df = ratio_rep_df[['id','return_ratio']]
    ratio_rep_dict = ratio_rep_df.set_index('id').T.to_dict()
    os.system('rm {i} {o}'.format(i=inside_bed, o=outside_bed))
    print('Temp files removed.\n')
    return(ratio_rep_dict) 
