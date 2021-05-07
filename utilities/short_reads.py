#import os, subprocess, sys

def star_mapping(index_dir, SR_fofn, output_dir):
    mapping_dir = output_dir + '/STAR_mapping'
    with open(SR_fofn) as fofn:
        for line in fofn:
            files = line.split(' ')
            sample_name = os.path.splitext(files[0])[-2].split('/')[-1]
            sample_prefix = mapping_dir + '/' + sample_name
            if not os.path.exists(sample_prefix + 'Log.final.out'): #I am replacing $index_name_STAR with output_dir + /STAR_index (not shure if this is the right argument to replace)
                print('Mapping for ', sample_name, ': in progress...')
                if len(files) == 1:
                    mapping = subprocess.Popen(['STAR', '--runThreadN',  '8', '--twopassMode', 'Basic', '--genomeDir', index_dir, '--readFilesIn', files[0], '--outFileNamePrefix', sample_prefix,'--alignSJoverhangMin', '8', '--alignSJDBoverhangMin', '1', '--outFilterType', 'BySJout', '--outSAMunmapped', 'Within', '--outFilterMultimapNmax', '20', '--outFilterMismatchNoverLmax', '0.04', '--outFilterMismatchNmax', '999', '--alignIntronMin', '20', '--alignIntronMax', '1000000', '--alignMatesGapMax', '1000000', '--sjdbScore', '1', '--genomeLoad', 'NoSharedMemory', '--outSAMtype', 'BAM', 'SortedByCoordinate'])
                    mapping.wait()
                    print('Mapping for ', sample_name, ': done.')
                else:
                    mapping = subprocess.Popen(['STAR', '--runThreadN', '8', '--twopassMode', 'Basic', '--genomeDir', index_dir, '--readFilesIn', files[0], files[1], '--outFileNamePrefix', sample_prefix,'--alignSJoverhangMin', '8', '--alignSJDBoverhangMin', '1', '--outFilterType', 'BySJout', '--outSAMunmapped', 'Within', '--outFilterMultimapNmax', '20', '--outFilterMismatchNoverLmax', '0.04', '--outFilterMismatchNmax', '999', '--alignIntronMin', '20', '--alignIntronMax', '1000000', '--alignMatesGapMax', '1000000', '--sjdbScore', '1', '--genomeLoad', 'NoSharedMemory', '--outSAMtype', 'BAM', 'SortedByCoordinate'])
                    mapping.wait()
                    print('Mapping for ', sample_name, ': done.')

def star(genome, SR_fofn, output_dir):
    fasta_genome = genome #should we check if it is fasta format?
    index_dir = output_dir + '/STAR_index'
    mapping_dir = output_dir + '/STAR_mapping'
    print('START')
    if not os.path.exists(mapping_dir):
        os.makedirs(mapping_dir)
    if not os.path.exists(index_dir):
        os.makedirs(index_dir)
        print('Running indexing...')
        indexing = subprocess.Popen(['STAR', '--runThreadN', '8', '--runMode', 'genomeGenerate', '--genomeDir', index_dir, '--genomeFastaFiles', fasta_genome], stdout=subprocess.PIPE)
        indexing.wait()
        print('Indexing done.')
    else:
        print('Index folder identified. Proceeding to mapping.')
    star_mapping(index_dir, SR_fofn, output_dir)
    return(mapping_dir)

#def kallisto_index ():
	
def kallisto_quantification(corrected_fasta, SR_fofn, output_dir):
    kallisto_output = output_dir + "/kallisto_output"
    if not os.path.exists(kallisto_output):
        os.makedirs(kallisto_output)
        kallisto_index = kallisto_output + '/kallisto_corrected_fasta.idx'
        run_kallisto = subprocess.Popen(['kallisto', 'index', '-i', kallisto_index, corrected_fasta])
        run_kallisto.wait()
    if not os.path.exists(kallisto_output + '/abundance.tsv'):
        os.makedirs(output_dir + '/STAR_mapping')
        index_dir=output_dir + '/STAR_index'
        star_mapping(index_dir, SR_fofn, output_dir)



#if args.SR_fofn and not args.coverage:
#	import utilities/short_reads as sr
#	print ('Short-read files provided.')
#	args.coverage = sr.star(args.genome, args.SR_fofn, args.dir)


#gtf = get_corr_filenames(args, dir=None)[0]???
def wobbler(gtf, mapping_dir): #mapping_dir: a directory where bam files are located
    replicates = []
    for filename in os.listdir(mapping_dir):
        if filename.endswith('sorted.bam'):
            replicates.append(filename)
    print('Replicates identified: '+str(replicates))
    os.system('grep -w "transcript" %s > temp.gtf' %gtf) 
    subprocess.call(["""samtools view -H %s | grep @SQ | sed 's/@SQ\tSN:\|LN://g' > chr_ordered.txt""" %(mapping_dir+'/'+replicates[0])], shell=True)
    subprocess.call(['cat chr_ordered.txt | cut -f1 > pattern.txt'], shell=True)
    subprocess.call(["""cat temp.gtf | awk '{if ($7 == "+") print $1 "\t" $4 "\t" $4+100 "\t" $12 "\t0\t"$7 ; else print $1 "\t" $5-100 "\t" $5 "\t" $12 "\t0\t"$7 }'>%s""" %'in_exons.bed'], shell=True)
    subprocess.call(["""cat temp.gtf | awk '{if ($7 == "+" && $4 < 101)  print $1 "\t1\t" $4 "\t" $12 "\t0\t"$7 ; else if ($7 == "+") print $1 "\t"$4-101"\t" $4 "\t" $12 "\t0\t"$7 ; else print $1 "\t" $5+1 "\t" $5+101 "\t" $12 "\t0\t"$7 }' > %s""" %'out_exons.bed'], shell=True)
    sorted_bed = """sed s/';'//g %s | sed s/'"'//g | sort -k1,1 -k2,2n -k3,3n > %s"""
    sorted_2_bed = 'while read -r line;  do grep $line -w %s >> %s; done < pattern.txt'
    bedtools_coverage= 'bedtools coverage -a %s -b %s -sorted -s -counts -f 0.1 -g chr_ordered.txt > %s'
    print('in/out_exons.bed files generated')
    subprocess.call([sorted_bed %('in_exons.bed', 'in_exons.sorted.bed')], shell=True)
    subprocess.call([sorted_bed %('out_exons.bed', 'out_exons.sorted.bed')], shell=True)
    subprocess.call([sorted_2_bed %('in_exons.sorted.bed', 'in_exons.sorted.2.bed')], shell=True)
    subprocess.call([sorted_2_bed %('out_exons.sorted.bed', 'out_exons.sorted.2.bed')], shell=True)
    print('in/out_exons.sorted.bed files generated')
    for i in ('in', 'out'):
        for k in range(len(replicates)):
            r=k+1
            bam_file=mapping_dir+'/'+replicates[k]
            subprocess.call([bedtools_coverage %(i+'_exons.sorted.2.bed', bam_file, "%sside_counts.rep%s.tsv"%(i, r))], shell=True)
            subprocess.call(["""cut -f7 %sside_counts.rep%s.tsv > %sside.rep%s.temp.counts"""%(i,r,i,r)], shell=True)
        subprocess.call(['paste *temp.counts > %sside.counts'%i], shell=True)
        subprocess.call(['cut -f4 %sside_counts.rep%s.tsv > replicate.counts'%(i,r)], shell=True)
        os.system('rm *temp.counts')
        col='$1'
        for k in range(len(replicates)):
            if (k+2)<=len(replicates):
                col = col+"+$%s"%str(k+2)
            else:
                break
        cmd = """cat %sside.counts | awk '{ print (%s)/%s}'>average.counts"""
        subprocess.call([cmd %(i, col, str(len(replicates)))], shell=True)
        subprocess.call(['paste replicate.counts average.counts > %sside_TSS.counts.tsv'%i], shell=True)
        os.system('rm %sside.counts average.counts'%i)
        print('%sside_TSS.counts.tsv generated'%i)


    #subprocess.call([counts %('in_exons.sorted.2.bed', bam, 'inside_TSS.counts.tsv')], shell=True)
    #subprocess.call([counts %('out_exons.sorted.2.bed', bam, 'outside_TSS.counts.tsv')], shell=True)
    #print('outside_TSS.counts.tsv files generated')
    #calculate ratios
    subprocess.call(["""paste inside_TSS.counts.tsv outside_TSS.counts.tsv | awk '{ print $1"\t"(($2+0.0001)/($4+0.0001)) }'>ratio_in_out_counts.tsv"""], shell=True)
    os.system('rm pattern.txt temp.gtf in_exons.bed in_exons.sorted.bed out_exons.bed out_exons.sorted.bed chr_ordered.txt in_exons.sorted.2.bed out_exons.sorted.2.bed inside_TSS.counts.tsv outside_TSS.counts.tsv')
    print('Temp files removed.\nRatios file "ratio_in_out_counts.tsv" generated!')

