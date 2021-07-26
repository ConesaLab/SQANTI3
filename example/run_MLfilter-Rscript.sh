# Command to generate ML filter outputs from SQANTI3 output
# Only includes mandatory arguments, all parameters set to default
# Run from SQANTI3 directory!

Rscript utilities/SQANTI3_MLfilter.R -c example/SQANTI3_output/UHR_chr22_classification.txt -o UHR_chr22 -d example/MLfilter_output 

