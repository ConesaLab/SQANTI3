# Command to run the ML filter using the SQANTI3 output                                                                                                                                                     
# Only includes mandatory arguments, all parameters set to default                                                                                                                                          
# Run from SQANTI3 directory!                                                                                                                                                                               

python sqanti3_filter.py ml --sqanti_class example/SQANTI3_QC_output/UHR_chr22_classification.txt -d example/MLfilter_output/ -o UHR_chr22 --filter_isoforms example/SQANTI3_QC_output/UHR_chr22_corrected.fasta 
