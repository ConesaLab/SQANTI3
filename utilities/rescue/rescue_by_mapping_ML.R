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
  optparse::make_option(c("-u", "--utilities_path"), type = "character",
                        help = "Full path to SQANTI3/utilities folder."),
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
                by = "mapping_hit") %>% 
      dplyr::rename(hit_POS_MLprob = "POS_MLprob")

    # add structural categories of candidates to mapping hits table
    mapping_hits <- mapping_hits %>% 
      dplyr::rename(isoform = "rescue_candidate") %>% 
      dplyr::left_join(classif %>% 
                         dplyr::select(isoform, structural_category), 
                       by = "isoform") %>% 
      dplyr::rename(rescue_candidate = "isoform") %>% 
      dplyr::rename(candidate_structural_category = "structural_category")
    
  
#### PERFORM RESCUE ####
    
    ## 1. Filter by ML probability threshold
    
      # filter mapping hits by probability threshold
      # and select hit based on max probability
      mapping_hits.max <- mapping_hits %>% 
        dplyr::filter(hit_POS_MLprob >= opt$threshold) %>% 
        dplyr::group_by(rescue_candidate) %>%
        dplyr::filter(hit_POS_MLprob == max(hit_POS_MLprob)) %>% 
        dplyr::ungroup()
      
    ## 2. Select only reference+unique rescued transcripts
      
      # select rescued transcripts from reference
      rescued_ref <- mapping_hits.max %>% 
        dplyr::filter(mapping_hit %in% probs.ref$isoform)
      
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
    
      # process automatic rescue result
      
          # find FSM rescued during automatic rescue to add to rescue table
          automatic_fsm <- classif %>% 
            dplyr::select(isoform, associated_transcript, 
                          structural_category) %>% 
            dplyr::right_join(automatic_ref_rescued %>% dplyr::rename(associated_transcript = "ref_transcript"), 
                              by = "associated_transcript")
          
          # add ML probabilities for automatic rescued references
          # and add SAM flag, rescue_result and exclusion_reason columns
          automatic_fsm <- automatic_fsm %>% 
            dplyr::left_join(probs.ref %>% dplyr::rename(associated_transcript = "isoform"), 
                                                          by = "associated_transcript") %>% 
            dplyr::mutate(sam_flag = NA, 
                          rescue_result = "rescued_automatic", 
                          exclusion_reason = NA)
          
          # reorder and rename
          automatic_fsm <- automatic_fsm %>% 
            dplyr::rename(rescue_candidate = "isoform", 
                          mapping_hit = "associated_transcript", 
                          hit_POS_MLprob = "POS_MLprob", 
                          candidate_structural_category = "structural_category") %>% 
            dplyr::relocate(sam_flag, .after = rescue_candidate) %>% 
            dplyr::relocate(hit_POS_MLprob, .after = mapping_hit)
          
      
      # include final rescue result in mapping hits table
      rescue_table <- mapping_hits %>% 
        dplyr::mutate(rescue_result = dplyr::if_else(
          mapping_hit %in% rescued_mapping_final$ref_transcript, 
          true = "rescued_mapping",
          false = "not_rescued"))
      
      # include exclusion reason for those not rescued
      rescue_table <- rescue_table %>% 
        dplyr::mutate(exclusion_reason = dplyr::case_when(
          # hits excluded due to insufficient probability, i.e. are not present in mapping_hits.max
          mapping_hit %in% mapping_hits.max$mapping_hit == FALSE ~ "ML_probability",
          
          # hits with high ML probability that constitute long read transcripts
          mapping_hit %in% mapping_hits.max$mapping_hit & 
            !(mapping_hit %in% probs.ref$isoform) ~ "long_read_transcript",
          
          # hits that are included in rescued_ref have both high probability and are not from long reads
          # and hits that are in isoform_assoc.tr have already been rescued or 
          # are represented by an LR transcript
          mapping_hit %in% rescued_ref$mapping_hit &
            mapping_hit %in% isoform_assoc.tr$associated_transcript ~ "reference_already_present"
        ))
      
      # join FSM/automatic rescue results
      rescue_table <- dplyr::bind_rows(rescue_table, automatic_fsm)
      
      # create best match column
      rescue_table <- rescue_table %>% 
        dplyr::group_by(rescue_candidate) %>% 
        dplyr::mutate(best_match_for_candidate = dplyr::case_when(
          # if there is a good matching reference transcript, set match column to ref
          any(exclusion_reason == "reference_already_present" |
                rescue_result == "rescued_mapping" |
                rescue_result == "rescued_automatic") ~ "reference_transcript",
          
          # if no good matching ref transcript was found but there is at least one 
          # good LR transcript, set match column to LR
          all(rescue_result == "not_rescued" & 
                exclusion_reason != "reference_already_present") &
            any(exclusion_reason == "long_read_transcript") ~ "long_read_transcript",
          
          # if none of the above is true, all hits were excluded due to ML probability
          # and match column is set to unknown (no match could be validated)
          all(exclusion_reason == "ML_probability") ~ "unknown"
        ))
      
      
      # create candidate - best match ID table
      
          # filter by ML probability (by rescue candidate groups)
          rescue_table.max <- rescue_table %>% 
            dplyr::filter(hit_POS_MLprob >= opt$threshold) %>% 
            dplyr::filter(hit_POS_MLprob == max(hit_POS_MLprob))
          
          # get rescue candidates with clear best match by probability
          match_unique.ids <- rescue_table.max %>% 
            dplyr::filter(dplyr::n() == 1) %>% 
            dplyr::select(rescue_candidate, mapping_hit) %>% 
            dplyr::rename(best_match_id = "mapping_hit")
            
          # find out which rescue candidates have ambiguity/ties
          rescue_ties <- rescue_table.max %>% 
            dplyr::summarize(matches = dplyr::n()) %>% 
            dplyr::filter(matches > 1) %>% 
            dplyr::select(rescue_candidate) %>% unlist
          
          # run
          rescue_table.ties <- rescue_table.max[
            rescue_table.max$rescue_candidate %in% rescue_ties,]
          
          # split table for candidates with primary alignments
          # and candidates with no primary alignments
          rescue_table.ties.prim <- rescue_table.ties %>%
            dplyr::filter(sam_flag == 0)
          
          # ids of candidate with no primary alignments
          nonprim <- setdiff(rescue_table.ties$rescue_candidate,
                             rescue_table.ties.prim$rescue_candidate)
          
          rescue_table.ties.nonprim <- rescue_table.ties[
            rescue_table.ties$rescue_candidate %in% nonprim,]
          
          # get best matches using primary alignments only
          match_tie.ids.prim <- rescue_table.ties.prim %>%
            dplyr::group_by(rescue_candidate) %>%
            dplyr::summarise(best_match_id = paste(mapping_hit, collapse=","))
          
          # get best matches using secondary alignments only
          match_tie.ids.nonprim <- rescue_table.ties.nonprim %>%
            dplyr::group_by(rescue_candidate) %>%
            dplyr::summarise(best_match_id = paste(mapping_hit, collapse=","))
          
          # merge matches
          match_tie.ids <- rbind(match_tie.ids.prim, match_tie.ids.nonprim)
          
          # join match tables
          match_ids <- dplyr::bind_rows(match_tie.ids, 
                                        match_unique.ids %>% dplyr::ungroup())
      
          
      # add match ID column to rescue table
      rescue_table <- rescue_table %>% 
        dplyr::left_join(match_ids, 
                         by = "rescue_candidate")
      
          # handle match ID col NAs caused by:
          # unknown best match cases
          rescue_table <- rescue_table %>% 
            dplyr::mutate(best_match_id = dplyr::if_else(best_match_for_candidate == "unknown", 
                                                         true = "unknown", 
                                                         false = best_match_id))
          # automatic rescue cases
          rescue_table <- rescue_table %>% 
            dplyr::mutate(best_match_id = dplyr::if_else(rescue_result == "rescued_automatic", 
                                                         true = mapping_hit,
                                                         false = best_match_id))
      
      # output rescue table
      readr::write_tsv(rescue_table,
                       file = paste0(opt$dir, "/", opt$output, 
                                     "_rescue_table.tsv"))
