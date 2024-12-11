#!/bin/bash Rscript
#
###################################################
#########     SQANTI3 AUTOMATIC RESCUE    #########
###################################################
#
# Author: Ángeles Arzalluz-Luque
# Contact: angeles.arzalluz@gmail.com
# Affiliation: Institute for Integrative Systems Biology, CSIC, Valencia, Spain
#
# DESCRIPTION:
# This script contains the code for the first module of the SQANTI3 rescue
# pipeline, which performs an automatic rescue of FSM-supported reference 
# transcripts that are lost during filtering. It outputs a list of long read
# transcripts that belong to other categories and are candidates for rescue,
# as well as list of reference transcripts that are ready to be rescued.
#


#### ARGUMENTS ####

# script argument list
option_list = list(
  optparse::make_option(c("-c","--sqanti_filter_classif"), type="character", default = NULL, 
                        help="SQANTI filter output classification file."),
  optparse::make_option(c("-o","--output"), type="character", default = "SQANTI3", 
                        help="Output file prefix."),
  optparse::make_option(c("-d","--dir"), type="character", 
                        help="Output directory."),
  optparse::make_option(c("-u", "--utilities_path"), type = "character",
                        help = "Full path to SQANTI3/utilities folder."),
  optparse::make_option(c("-g", "--refGTF"), type = "character",
                        help = "Full path to reference transcriptome GTF used when
                        running SQANTI3."), 
  optparse::make_option(c("-e", "--rescue_mono_exonic"), type = "character", default = "all",
                        help = "Whether or not to include mono-exonic artifacts 
                        in the rescue. Options include: 'none', 'fsm' and 'all' (default)."),
  optparse::make_option(c("-m", "--mode"), type = "character", default = "automatic",
                        help = "If 'automatic' (default), only automatic rescue of FSM 
                        artifacts will be performed. If 'full', rescue will include mapping 
                        of ISM, NNC and NIC artifacts to find potential replacement isoforms.")
  )


# Parse and handle provided arguments
opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser) # list of the args


#### PREPARATION ####

    # import pipe operator
    require(magrittr)

    # message to indicate which mode is activated
    message("\n---------------------------------------------------------------")
    message("\n\t\tINITIATING SQANTI3 RESCUE...\n")
    message("\n---------------------------------------------------------------")
    
    if(opt$mode == "automatic"){
      message("\n\t--mode automatic:")
      message("\n\t\tAutomatic rescue mode selected (default).\n") 
      message("\n\t\tRescue will be performed only for artifact FSM transcripts.\n")
    } else{
      message("\n\t--mode full:")
      message("\n\t\tFull rescue mode selected!\n") 
      message("\n\t\tAutomatic rescue activated for artifact FSM transcripts.")
      message("\n\t\tAdditional rescue steps will be performed for ISM, NIC and NNC artifacts.\n")
    }
   
    # message for loading classification file
    message("\n---------------------------------------------------------------")
    message("\n\tREADING FILTER CLASSIFICATION FILE...\n")

    # read in classification
    classif <- readr::read_tsv(opt$sqanti_filter_classif)
    
    message("\n---------------------------------------------------------------")
    

####----------------- AUTOMATIC RESCUE OF FSM ----------------####

message("\n---------------------------------------------------------------")
message("\n\tPERFORMING AUTOMATIC RESCUE...\n")
message("\n---------------------------------------------------------------")

# print warning message regarding -e option
message(paste0("\n\t***NOTE: you have set -e ", opt$rescue_mono_exonic, ":"))
if(opt$rescue_mono_exonic == "all"){
  message("\n\t\tAll mono-exonic artifact transcripts will be considered for rescue.")
}else if(opt$rescue_mono_exonic == "fsm"){
  message("\n\t\tMono-exonic artifact transcripts will only be considered for rescue if they are FSM.")
}else{
  message("\n\t\tAll mono-exonic artifact transcripts will be excluded from the rescue.")
}

## MONO-EXONS ##
# rescue mono-exon FSM transcripts (if indicated)

if(opt$rescue_mono_exonic %in% c("all", "fsm")){

    message("\n\tRescuing references associated to mono-exon FSM...")

    # filter classification
    classif_mono <- classif %>% 
      dplyr::filter(structural_category == "full-splice_match" & 
                      exons == 1 & filter_result == "Artifact")
    
    # select associated transcripts
    rescue_mono <- classif_mono %>% 
      dplyr::select(associated_transcript) %>% 
      unique %>% 
      dplyr::rename(isoform = "associated_transcript")
}


## FSM and ISM ##
# analyze lost reference transcripts and 
# see whether they are supported by ISM/FSM
    
    # filter classification
    classif_ism_fsm <- classif %>% 
      dplyr::filter(structural_category %in% c("full-splice_match", 
                                               "incomplete-splice_match") &
                      exons > 1)
    
    # add ISM mono-exons if indicated
    if(opt$rescue_mono_exonic == "all" & 
       opt$mode  == "full"){
      
      message("\n\tIncluding mono-exon ISM as rescue candidates...")
      
      classif_ism_mono <- classif %>% 
        dplyr::filter(structural_category == "incomplete-splice_match" &
                        exons == 1)
      classif_ism_fsm <- dplyr::bind_rows(classif_ism_fsm,
                                          classif_ism_mono)
    }
    
    # find all reference IDs in associated_transcript column
    all_ref <- classif_ism_fsm %>% 
      dplyr::select(associated_transcript) %>% unique %>% unlist
    
    # check reference IDs not represented by isoforms (lost in filtering)
    isoform_ref <- classif_ism_fsm %>% 
      dplyr::filter(filter_result == "Isoform") %>% 
      dplyr::select(associated_transcript) %>% unique %>% unlist
    
    lost_ref <- all_ref[!(all_ref %in% isoform_ref)]
    
    
    # perform rescue for lost reference transcripts
    message("\n\tFinding FSM-supported reference transcripts lost after filtering...")
    
    source(paste0(opt$utilities_path, "/rescue/rescue_aux_functions.R"))
    rescue <- purrr::map_df(lost_ref, rescue_lost_reference, classif_ism_fsm)

    # separate result into reference transcripts and ISM
    rescue_ism <- rescue %>% 
      dplyr::filter(isoform %in% classif_ism_fsm$isoform)
    
    rescue_ref <- rescue %>% 
      dplyr::filter(isoform %in% classif_ism_fsm$associated_transcript)
    
    
    # write out reference transcripts that are automatically rescued
    # add mono-exons if indicated
    if(opt$rescue_mono_exonic %in% c("all", "fsm")){
      rescue_auto <- dplyr::bind_rows(rescue_mono, rescue_ref)
    }else{
      rescue_auto <- rescue_ref
    }
    
    readr::write_tsv(rescue_auto, 
                     col_names = FALSE,
                     file = paste0(opt$dir, "/", opt$output, 
                                   "_automatic_rescued_list.tsv"))
    
    # write out rescue table when mode == "automatic"
    # includes all artifacts (ID/category) considered during automatic rescue
    # and the IDs of the rescued reference transcripts
    if(opt$mode == "automatic"){
      
      # get ISM/FSM transcripts associated to all automatically rescued references
      rescue_table <- dplyr::filter(classif_ism_fsm, 
                                    associated_transcript %in% rescue_auto$isoform) %>% 
        dplyr::select(isoform, associated_transcript, structural_category) %>% 
        dplyr::rename(artifact = "isoform", 
                      rescued_transcript = "associated_transcript")
      
      # write table
      readr::write_tsv(rescue_table,
                       file = paste0(opt$dir, "/", opt$output,
                                     "_automatic_rescue_table.tsv"))
    }
    
    
  message("\n\tAutomatic rescue finished successfully!")
  message("\n\tReference transcripts output by automatic rescue were written to output file:")
  message(paste0("\n\t\t", opt$dir, "/", opt$output, 
                 "_automatic_rescued_list.tsv"))
  message(paste0("\n\t\tTotal transcripts rescued from reference : ", 
                 rescue_auto %>% nrow))
    if(opt$rescue_mono_exonic %in% c("all", "fsm")){
      message(paste0("\n\t\t - From mono-exon FSM: ", rescue_mono %>% nrow))
    }else{
      message("\n\t\t - From mono-exon FSM: 0")
    }
  message(paste0("\n\t\t - From FSM artifacts: ", 
                 rescue_ref %>% nrow))

message("\n---------------------------------------------------------------")

    

####----------------- RESCUE CANDIDATES ----------------####

## NNC and NIC ###
# get NIC and NNC artifacts and join ISM not associated with rescued references
# only when mode == "full" (!!)

if(opt$mode == "full"){
  
  message("\n---------------------------------------------------------------")
  message("\n\tFINDING RESCUE CANDIDATES...\n")
  message("\n---------------------------------------------------------------")
  message("\n\tRescue candidates: artifact transcripts to be used for rescue.\n")
  
  # filter classification to get NIC and NNC artifacts
  rescue_novel <- classif %>% 
    dplyr::filter(filter_result == "Artifact" &
                    structural_category %in% c("novel_in_catalog", 
                                               "novel_not_in_catalog"))
  
  # exclude mono-exonic if indicated
  if(opt$rescue_mono_exonic %in% c("fsm", "none")){
    
    all_novel <- nrow(rescue_novel)
    
    rescue_novel <- rescue_novel %>% 
      dplyr::filter(exons > 1)
    
    message(paste0("\n\tExcluding ", 
                   all_novel - nrow(rescue_novel),
                   " mono-exonic novel transcripts from rescue candidate list."))
  }
  
  # write out rescue candidates (novel and ISM)
  rescue_candidates <- dplyr::bind_rows(rescue_ism, 
                                        rescue_novel %>% dplyr::select(isoform))
  
  # remove mono-exons if instructed
  readr::write_tsv(rescue_candidates, 
                   col_names = FALSE,
                   file = paste0(opt$dir, "/", opt$output,
                                 "_rescue_candidates.tsv"))
  
  # print messages and summary
  message("\n\tRescue candidate search finished successfully!")
  message("\n\tArtifact transcripts flagged as rescue candidates were written to output file:")
  message(paste0("\n\t\t", opt$dir, "/", opt$output, "_rescue_candidates.tsv"))
  message("\n\tRescue candidate summary")
  message(paste0("\n\t\t ISM: ", rescue_ism %>% nrow))
  message(paste0("\n\t\t NIC: ", rescue_novel %>% 
                   dplyr::filter(structural_category == "novel_in_catalog") %>% 
                   nrow))
  message(paste0("\n\t\t NNC: ", rescue_novel %>% 
                   dplyr::filter(structural_category == "novel_not_in_catalog") %>% 
                   nrow))
  message(paste0("\n\t Total: ", rescue_candidates %>% nrow))
  
  message("\n---------------------------------------------------------------")
  
  
}
# END OF MODE == "FULL" CONDITION FOR NIC/NNC CANDIDATE SEARCH #


    
####----------------- RESCUE TARGETS ----------------####

# Isoforms passing filter from long-read transcriptome and reference transcripts
# only when mode == "full" (!!)

if(opt$mode == "full"){
  
  message("\n---------------------------------------------------------------")
  message("\n\tRETRIEVING RESCUE TARGETS...\n")
  message("\n---------------------------------------------------------------")
  message("\n\t Rescue targets: validated LR or reference isoforms that could replace an artifact from the same gene.\n")
  
  # obtain rescue targets for defined rescue candidates
  
  # get all genes with rescue candidates
  message("\n\t Retrieving target genes...\n")
  
  rescue_candidates <- rescue_candidates %>% 
    dplyr::left_join(classif %>% dplyr::select(isoform, associated_gene), 
                     by = "isoform")
  
  target_genes <- rescue_candidates$associated_gene %>% unique
  
  # get target isoforms for rescue from long read transcriptome
  message("\n\t Finding target isoforms from long read transcriptome...\n")
  
  rescue_targets.LR <- classif %>% 
    dplyr::filter(filter_result == "Isoform" &
                    associated_gene %in% target_genes) %>% 
    dplyr::select(isoform)
  
  # get targets from reference transcriptome
  message("\n\t Finding target isoforms from reference transcriptome...\n")
  
  reference_ids <- BUSpaRse::tr2g_gtf(opt$refGTF, 
                                      gene_version = NULL,
                                      get_transcriptome = FALSE, 
                                      out_path = opt$dir, 
                                      save_filtered_gtf = FALSE,
                                      write_tr2g = FALSE,
                                      chrs_only = FALSE)
  
  rescue_targets.ref <- reference_ids %>% 
    dplyr::filter(gene %in% target_genes) %>% 
    dplyr::select(transcript) %>% 
    dplyr::rename(isoform = "transcript")
  
  # join both target lists
  rescue_targets <- dplyr::bind_rows(rescue_targets.LR,
                                     rescue_targets.ref)
  
  # write to file
  readr::write_tsv(rescue_targets, 
                   col_names = FALSE,
                   file = paste0(opt$dir, "/", opt$output,
                                 "_rescue_targets.tsv"))
  
  # print messages and summary
  message("\n\tRescue target search finished successfully!")
  message("\n\tValidated isoforms to be used as rescue targets were written to output file:")
  message(paste0("\n\t\t", opt$dir, "/", opt$output, "_rescue_targets.tsv"))
  message("\n\tRescue target summary")
  message(paste0("\n\t\t LR transcriptome: ", rescue_targets.LR %>% nrow))
  message(paste0("\n\t\t Reference transcriptome: ", rescue_targets.ref %>% nrow))
  message(paste0("\n\t Total: ", rescue_targets %>% nrow))
  message("\n---------------------------------------------------------------")
  
}

# END OF MODE == "FULL" CONDITION FOR TARGET SEARCH #


