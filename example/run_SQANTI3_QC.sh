## This is the command used to create the SQANTI3 output for the example data

python sqanti3_qc.py example/UHR_chr22.gtf example/gencode.v38.basic_chr22.gtf example/GRCh38.p13_chr22.fasta --CAGE_peak data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed --polyA_motif_list data/polyA_motifs/mouse_and_human.polyA_motif.txt -o UHR_chr22 -d example/SQANTI3_QC_output -fl example/UHR_abundance.tsv --short_reads example/UHR_chr22_short_reads.fofn --cpus 4 --report both
