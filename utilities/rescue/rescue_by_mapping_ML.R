#!/bin/bash Rscript
#
####################################################
#########     SQANTI3 RESCUE BY MAPPING    #########
####################################################
#
# Author: Ángeles Arzalluz-Luque
# Contact: angeles.arzalluz@gmail.com
# Affiliation: Institute for Integrative Systems Biology, CSIC, Valencia, Spain
#
# DESCRIPTION:
# This script contains the code for the third module of the SQANTI3 rescue
# pipeline, which integrates rescue candidate alignment results (to their same-gene
# counterparts, both from the reference and from the long read defined transcriptome). 
# Using the ML probabilities obtained for all rescue targets (both reference 
# and long read-defined transcripts), it selects the best rescue target from 
# those that were retrieved during alignment.
#

#### ARGUMENTS ####

# script argument list
option_list = list(
  optparse::make_option(c("-c","--sqanti_MLclassif"), type = "character", default = NULL, 
                        help = "SQANTI ML output classification file."),
  optparse::make_option(c("-o","--output"), type = "character", default = "SQANTI3", 
                        help = "Output file prefix."),
  optparse::make_option(c("-d","--dir"), type = "character", 
                        help="Output directory."),
  optparse::make_option(c("-m", "--mapping_hits"), type = "character",
                        help = "Path to file containing artifact isoform pairs 
                        (rescue candidates and targets) obtained during alignment."),
  optparse::make_option(c("-r", "--reference_MLprob"), type = "character",
                        help = "Path to file containing reference transcriptome
                        ML probabilities (obtained after running trained random 
                        forest classifier on reference transcriptome)."),
  optparse::make_option(c("-j","--threshold"), type = "numeric", default = 0.7,
                        help = "Default: 0.7. Machine learning probability 
                        threshold to filter elegible rescue targets (mapping hits).")
)


# Parse and handle provided arguments
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser) # list of the args
opt$threshold <- as.numeric(opt$threshold)


#### PREPARATION ####

    # import pipe operator
    require(magrittr)
  
    # load mapping hits obtained from SAM
    mapping_hits <- readr::read_tsv(opt$mapping_hits, 
                                    col_names = c("rescue_candidate", 
                                                  "sam_flag", 
                                                  "mapping_hit"))
    
    # load ML filter probabilities obtained for reference transcripts
    probs.ref <- readr::read_tsv(opt$reference_MLprob)
    
    # load ML probabilities obtained for long read transcriptome
    # in previous run of the ML filter
    classif <- readr::read_tsv(opt$sqanti_MLclassif)
    
    probs.LR <- classif %>% 
      dplyr::select(isoform, POS_MLprob) %>% 
      dplyr::filter(!is.na(POS_MLprob))
    
    # join both reference and LR ML probabilities
    probs <- dplyr::bind_rows(probs.ref, probs.LR)
    
    # add ML probabilities of mapping hits to mapping hits table
    mapping_hits <- mapping_hits %>% 
      dplyr::left_join(probs %>% 
                         dplyr::rename(mapping_hit = "isoform"), 
                by = "mapping_hit")

    # add structural categories of candidates to mapping hits table
    mapping_hits <- mapping_hits %>% 
      dplyr::rename(isoform = "rescue_candidate") %>% 
      dplyr::left_join(classif %>% 
                         dplyr::select(isoform, structural_category), 
                       by = "isoform") %>% 
      dplyr::rename(rescue_candidate = "isoform")
    
  
#### PERFORM RESCUE ####
    
    ## 1. Filter by ML probability threshold
    
      # filter mapping hits by probability threshold
      # and select hit based on max probability
      mapping_hits.max <- mapping_hits %>% 
        dplyr::filter(POS_MLprob >= opt$threshold) %>% 
        dplyr::group_by(rescue_candidate) %>%
        dplyr::filter(POS_MLprob == max(POS_MLprob)) %>% 
        dplyr::ungroup()
      
    ## 2. Select only reference+unique rescued transcripts
      
      # select rescued transcripts from reference
      rescued_ref <- mapping_hits.max %>% 
        dplyr::filter(stringr::str_detect(mapping_hit, 
                                          "PB.", negate = TRUE))
      
    ## 3. Remove reference transcripts already represented in transcriptome
    ##    to avoid introducing redundancy as a result of the rescue
      
      # retrieve all reference transcripts (associated_transcript)
      # that are already represented by an isoform
      isoform_assoc.tr <- classif %>%
        dplyr::filter(filter_result == "Isoform" & 
                        associated_transcript != "novel") %>% 
        dplyr::select(associated_transcript)
      
      # include those that were retrieved in automatic rescue
      automatic_ref_rescued <- readr::read_tsv(paste0(opt$dir, "/", opt$output, 
                                                      "_automatic_rescued_list.tsv"),
                                               col_names = "associated_transcript")
      
      isoform_assoc.tr <- dplyr::bind_rows(isoform_assoc.tr, 
                                           automatic_ref_rescued) %>% unique
      
      # find truly rescued references (not represented by any isoform)
      rescued_mapping_final <- rescued_ref %>% 
        dplyr::filter(!(mapping_hit %in% 
                        isoform_assoc.tr$associated_transcript)) %>% 
        dplyr::select(mapping_hit) %>% 
        dplyr::rename(ref_transcript = "mapping_hit") %>% 
        unique()
      
      # make compatible colnames
      automatic_ref_rescued <- automatic_ref_rescued %>% 
        dplyr::rename(ref_transcript = "associated_transcript")
      
      # generate final list of rescued transcripts
      rescued_final <- dplyr::bind_rows(automatic_ref_rescued,
                                        rescued_mapping_final)
      
      
#### WRITE OUTPUTS ####
      
      # output rescue inclusion list
      readr::write_tsv(rescued_final, 
                       col_names = FALSE,
                       file = paste0(opt$dir, "/", opt$output, 
                                     "_rescue_inclusion-list.tsv"))
    
      # include final rescue result in mapping hits table
      mapping_hits <- mapping_hits %>% 
        dplyr::mutate(rescue_result = dplyr::case_when(
          mapping_hit %in% automatic_ref_rescued$ref_transcript ~ "rescued_automatic",
          mapping_hit %in% rescued_mapping_final$ref_transcript ~ "rescued_mapping",
          mapping_hit %in% rescued_final$ref_transcript == FALSE ~ "not_rescued"),
        exclusion_reason = dplyr::case_when(
          mapping_hit %in% rescued_final$ref_transcript ~ NA,
          mapping_hit %in% mapping_hits.max$mapping_hit == FALSE ~ "MLprob",
          mapping_hit %in% mapping_hits.max$mapping_hit & 
            str_detect(mapping_hit, "PB.") ~ "LR",
          mapping_hit %in% rescued_ref$mapping_hit &
            mapping_hit %in% isoform_assoc.tr$associated_transcript ~ "reference_already_present"
        ))
      
      # output rescue table
      readr::write_tsv(mapping_hits,
                       file = paste0(opt$dir, "/", opt$output, 
                                     "_rescue_table.tsv"))
