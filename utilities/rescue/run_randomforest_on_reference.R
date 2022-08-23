#!/bin/bash Rscript
#
#######################################################################
#########     SQANTI3 REFERENCE ML CLASSIFIER RUN (RESCUE)    #########
#######################################################################
#
# Author: Ángeles Arzalluz-Luque
# Contact: angeles.arzalluz@gmail.com
# Affiliation: Institute for Integrative Systems Biology, CSIC, Valencia, Spain
#
# DESCRIPTION:
# This script contains the code for re-running a pre-trained ML filter classifier
# on the reference transcriptome used for a SQANTI3 run. To achieve this, SQANTI3
# QC is previously run on said reference transcriptome using the same orthogonal
# data as in the real SQANTI3 QC run. The script outputs a two-column table containing
# the transcript IDs for reference isoforms that were selected as rescue targets
# and the positive ML probability obtained as a result of running the classifier
# (POS_MLprob, i.e. the probability that the transcript is truly an isoform based
# on the attributes previously learned by the random forest classifier).
#


#### ARGUMENTS ####

# script argument list
option_list = list(
  optparse::make_option(c("-c","--classification"), type="character", default = NULL, 
                        help="Classification file obtained as a result of 
                        running SQANTI3 QC on the reference transcriptome."),
  optparse::make_option(c("-o","--output"), type="character", default = "SQANTI3", 
                        help="Output file prefix."),
  optparse::make_option(c("-d","--dir"), type="character", 
                        help="Output directory."),
  optparse::make_option(c("-r", "--randomforest"), type = "character",
                        help = "Full path to the randomforest.RData object
                        obtained during previous run of SQANTI3 ML filter.")
)

# Parse and handle provided arguments
opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser) # list of the args


message("\n-------------------------------------------------------------------------")
message("\n\tRUN PRE-TRAINED SQ3 ML CLASSIFIER RUN ON REFERENCE TRANSCRIPTOME\n")
message("\n-------------------------------------------------------------------------")


#### DATA PREPARATION ####

message("\n\tReading SQANTI3 classification file...\n")

    # load data (exactly as in ML filter)
    classification <- read.table(opt$classification, 
                                 sep ="\t", header= TRUE, as.is =TRUE)
    
    # set rownames
    classification <- as.data.frame(classification)
    rownames(classification) = classification$isoform
    classification = classification[,-which(colnames(classification) == "isoform")]


  
#### DATA CLEANING FOR ML ####
## Recreate all column-specific modifications used in ML filter to allow
## running the random forest algorithm on the SQ3 classification file
        
message("\n\tPerforming ML filter data cleaning on classification...\n")

    ## Replace NAs
    
    # Columns with NAs
    NA_columns <- c("within_CAGE_peak", 'n_indels', "n_indels_junc", "FL",
                    "predicted_NMD", "min_sample_cov", "min_cov", "ratio_exp", "bite", 
                    "diff_to_gene_TSS", "diff_to_gene_TTS" , "dist_to_polyA_site", 
                    "dist_to_CAGE_peak", 'within_polyA_site', "polyA_dist")
    
    replacement.na <- c(0, 0, 0, 0, "non_coding",0, 0,0, FALSE, 
                        -11000, -11000, -11000, -11000, FALSE, -11000)
    
    for (i in 1:length(NA_columns)) {
      sel.column <- which(colnames (classification) == NA_columns [i])
      classification[which(is.na(classification[,sel.column])), sel.column] <- replacement.na[i]
    }
    
    # Special case for sd_cov: If all NAs, replace by 0. If there are NA values, 
    # replace with  median value of sd_cov in the dataset
    if (all(is.na(classification$sd_cov))){
      classification$sd_cov = 2
    } else{
      medt2 = median(as.numeric(classification[!is.na(classification$sd_cov), "sd_cov"]))
      classification[is.na(classification$sd_cov),"sd_cov"] <- medt2
    }
    
    
    ## Handle column types
    
    # Convert in factors columns with categorical variables
    categorical <- c("FSM_class", "coding", "bite", "within_CAGE_peak", 
                     "polyA_motif_found" , "within_polyA_site", "predicted_NMD")
    for (x in categorical){
      classification[,x] <- as.factor(classification[,x])
    }
    
    # Convert in integers columns with numerical variables
    integers <- c("diff_to_gene_TSS", "diff_to_gene_TTS", "min_sample_cov",
                  "min_cov", "ratio_exp" , "polyA_dist", "dist_to_CAGE_peak",
                  "FL", "length", "exons")
    for (x in integers){
      classification[,x] <- as.integer(classification[,x])
    }
    

  
#### ENSURE COMPATIBILITY WITH PREVIOUS CLASSIFIER ####
## The trainingData slot stores the training dataframe as supplied to the random
## forest algorithm. New data must have exactly the same columns (name and format)
## that the previous classifier was trained under.
    
message("\n\tChecking column compatibility with pre-trained random forest classifier...\n")

    require(magrittr)
  
    ## Load pre-trained randomforest.RData object generated in ML filter run
    randomforest <- readRDS(opt$randomforest)
    
    ## Filter columns that were not used during training
    model_cols <- randomforest$trainingData %>% colnames
    classification <- classification[, model_cols[-length(model_cols)]]
    
    ## Check for NAs
    message("\n\tValidating columns used in prediction...")
    detect_na <- purrr::map_lgl(classification, ~(is.na(.) %>% any))
    message("\n\tColumn-level NA check:")
    print(detect_na)
    
    ## Check column types
    col_types <- purrr::map_chr(classification, class)
    message("\n\tColumn type check:")
    print(col_types)
    
    
#### RUN ML ####
    
message("\n\tRunning random forest classifier on reference transcriptome...\n")
    # trick to equilize classes of classification and training set
    classification = rbind(randomforest$trainingData[1,-dim(randomforest$trainingData)[2]] , classification)
    classification = classification[-1,]
    # apply classifier to data
    isoform_predict = predict(randomforest, classification, type = 'prob')
    
    rownames(isoform_predict) <- rownames(classification)
    isoform_predict <- isoform_predict %>% 
      tibble::rownames_to_column("isoform") %>% 
      tibble::as_tibble() %>% 
      dplyr::select(isoform, POS) %>% 
      dplyr::rename(POS_MLprob = "POS")
    
    readr::write_tsv(isoform_predict, 
                     file = paste0(opt$dir, "/", opt$output, 
                                   "_reference_isoform_predict.tsv"))

message("\n\tRandom forest classifier run successfully!\n")
message("\n\tResulting isoform/artifact probabilities written to output file:\n")
message(paste0("\n\t\t", opt$dir, "/", opt$output, 
               "_reference_isoform_predict.tsv", "\n"))
message("\n-------------------------------------------------------------------------")

