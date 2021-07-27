#!/bin/bash Rscript
#
#-------------------------------------------------------------------------------
#                       _________________________
#
#                         SQANTI3 FILTER REPORT
#                       _________________________
#
# Author: Ángeles Arzalluz-Luque
# 
# Contact: angeles.arzalluz@gmail.com
#
# Affiliation: Institute for Integrative Systems Biology, CSIC, Valencia, Spain
#
# Last updated: 13/July/2021
#
#-------------------------------------------------------------------------------



####------------------------- FILTER REPORT ARG LIST -----------------------####

#### Define script arguments ####

option_list <- list(
  optparse::make_option(c("-d","--dir"), type = "character", 
                        help = "Output/input directory - must be same as SQ3 filter."),
  optparse::make_option(c("-o","--output"), type = "character", default = "SQANTI3", 
                        help = "Output report file prefix."),
  optparse::make_option(c("f", "--filter_type"), type = "character",
                        help = "Type of SQ3 filter that was run to generate input.
                        Must be one of the following: 'ml' or 'rules'.")
)

# Parse arguments
opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser)




####--------------------------- INPUTS & THEME -----------------------------####

#### Initialize script and load input data ####

message("\n-------------------------------------------------")
message("\n \t SQANTI3 Machine Learning filter")
message("\n--------------------------------------------------")

# Import pipe operator
require(magrittr)


# Read files in output/input directory
filter_outfiles <- dir(opt$dir)
paths <- paste0(opt$dir, "/", filter_outfiles)


# Detect files depending on filter type

if(opt$filter_type == "ml"){
  
  # Detect path and load ML output classification
  message("\nReading ML result classification table...")
  
  which_classif <- stringr::str_detect(paths, "MLresult_classification")
  path_classif <- paths[which_classif]
  classif <- readr::read_tsv(path_classif, 
                             col_types = readr::cols(exons = readr::col_integer(),
                                                     ref_exons = readr::col_integer()))
  
  
  # Detect path and load variable importance table
  message("\nReading classifier variable importance table...")
  
  which_imp <- stringr::str_detect(paths, "variable-importance_table")
  path_imp <- paths[which_imp]
  imp <- readr::read_tsv(path_imp, col_names = c("variable", "importance"))
  
  # format imp as factor
  imp <- imp %>% dplyr::mutate(variable = factor(variable) %>% 
                                 forcats::fct_reorder(importance))
  
}


# Format category factor for classification file
classif <- classif %>% 
  dplyr::mutate(structural_category = factor(structural_category) %>% 
                  forcats::fct_infreq() %>% 
                  forcats::fct_recode(ISM = "incomplete-splice_match",
                                      FSM = "full-splice_match",
                                      NNC = "novel_not_in_catalog",
                                      NIC = "novel_in_catalog",
                                      Intergenic = "intergenic",
                                      Antisense = "antisense",
                                      Genic = "genic",
                                      Fusion = "fusion",
                                      `Genic \nintron` = "genic_intron"))


#### Define plot theme ####

# Load ggplot2
require(ggplot2)

# Set theme parameters (from SQANTI3_report.R)
sq_theme <- theme_classic(base_family = "Helvetica") +
  theme(plot.title = element_text(lineheight=.4, size=15, hjust = 0.5)) +
  theme(plot.margin = unit(c(2.5,1,1,1), "cm")) +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=13),
        axis.text.x  = element_text(size=12),
        axis.title.y = element_text(size=13),
        axis.text.y  = element_text(vjust=0.5, size=12) ) +
  theme(legend.text = element_text(size = 11), 
        legend.title = element_text(size=12), 
        legend.key.size = unit(0.5, "cm"))

theme_set(sq_theme)


# Set category palette
cat.palette = c(FSM="#6BAED6", ISM="#FC8D59", NIC="#78C679", 
                NNC="#EE6A50", Genic="#969696", Antisense="#66C2A4", 
                Fusion="goldenrod1", Intergenic = "darksalmon", `Genic \nintron`="#41B6C4")




####-------------------------- FILTER REPORT PLOTS ------------------------####

#### Common filter plots ####
    
# Isoforms and artifacts by category
category_summary <- classif %>% 
  dplyr::group_by(structural_category, filter_result) %>% 
  dplyr::summarize(n = dplyr::n()) %>%
  dplyr::mutate(percent = n/sum(n))
    
      # Plot totals
      cat_totals <- ggplot(category_summary, 
                            aes(x = structural_category, y = n)) + 
        geom_bar(aes(alpha = filter_result, fill = structural_category), 
                  stat = "identity",
                  width = 0.8, color = "black", position = "dodge") +
        labs(x = "Structural category", y = "Transcript no.") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        scale_fill_manual("", values = cat.palette, guide = "none") +
        scale_alpha_manual("Filter result", values = c(0.3, 1))
        
      # Plot percentages
      cat_percent <- ggplot(category_summary, 
                            aes(x = structural_category, y = percent)) +
        geom_bar(aes(alpha = filter_result, fill = structural_category), 
                  stat = "identity",
                  width = 0.8, color = "black", position = "stack") +
        labs(x = "Structural category", y = "Transcript %") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        scale_y_continuous(labels = scales::percent_format()) +
        scale_fill_manual("", values = cat.palette, guide = "none") +
        scale_alpha_manual("Filter result", values = c(0.3, 1))

#### END common filter plots ####
      
      
      
#### ML filter plots ####
if(opt$filter_type == "ml"){
    
    ## Reason for calling artifacts: ML, intra-priming or both
    classif_artifacts <- classif %>% 
      dplyr::filter(filter_result == "Artifact") %>% 
      dplyr::mutate(reason = dplyr::case_when(
        intra_priming == TRUE & ML_classifier == "Negative" ~ "Both",
        intra_priming == TRUE | (ML_classifier == "Positive" | is.na(ML_classifier)) ~ "Intra-priming",
        intra_priming == FALSE & ML_classifier == "Negative" ~ "ML"))
    
    artifact_summary <- classif_artifacts %>% 
      dplyr::group_by(structural_category, reason) %>% 
      dplyr::summarize(n = dplyr::n()) %>% 
      dplyr::mutate(percent = n/sum(n))
    
        # totals
        artifact_totals <- ggplot(artifact_summary) +
          ggtitle("Reason to flag transcripts as artifacts, by category") +
          geom_bar(aes(x = structural_category, y = n, fill = reason), 
                   stat = "identity", position = "dodge", 
                   color = "black", width = 0.8) +
          labs(x = "Structural category", y = "No. of transcripts") +
          theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
          scale_fill_manual("Reason", values = c("springgreen4", "darkgoldenrod2", 
                                                 "steelblue3"))
        
        # percentage
        artifact_percent <- ggplot(artifact_summary) +
          ggtitle("Reason to flag transcripts as artifacts, by category") +
          geom_bar(aes(x = structural_category, y = percent, fill = reason), 
                   stat = "identity", position = "stack", 
                   color = "black", width = 0.8) +
          labs(x = "Structural category", y = "% Transcripts") +
          theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
          scale_fill_manual("Reason", values = c("springgreen4", "goldenrod2", 
                                                 "palevioletred3")) +
          scale_y_continuous(labels = scales::percent_format())
          
    
    ## Variable importance in classifier
    var_imp <- ggplot(imp) +
      ggtitle("Variable importance in Random Forest classifier") +
      geom_bar(aes(x = variable, y = importance), stat = "identity",
               width = 0.7, color = "black", fill = "cadetblue3") +
      labs(x = "SQANTI3 variables", y = "Importance") +
      coord_flip()
    
    
    ## Variables used in ML: values for isoforms and artifacts by category
    source(paste0(getwd(), "/utilities/compare_MLvariables.R"))
    
    var_compare <- purrr::map2(imp$variable, imp$importance,
                               ~compare_MLvariables(classif, .x, .y))
    
}
#### END ML filter plots ####

      
      
#### Intra-priming plots ####

# Intra-priming by category and exon number
ip_summary <- classif %>% 
    dplyr::select(structural_category, exons, intra_priming) %>% 
    dplyr::mutate(exon_classif = dplyr::if_else(exons > 1, 
                                                true = "Multi-exon", 
                                                false = "Mono-exon"))
      
ip_byexons <- ip_summary %>% 
    dplyr::group_by(structural_category, exon_classif, intra_priming) %>% 
    dplyr::summarize(n = dplyr::n()) %>% 
    dplyr::mutate(percent = n/sum(n))
  
  
    # A % summary
      ggplot(classif) +
        ggtitle("Percentage of A's downstream of TTS") +
        geom_density(aes(x = perc_A_downstream_TTS, fill = intra_priming), 
                     alpha = 0.5) +
        scale_fill_manual("Intra-priming prediction", 
                          values = c("violetred3", "gray60"))
      
      
    # Intra-priming by structural_category and exon_classif
        # totals
        ip_totals <- ggplot(ip_summary) +
          geom_bar(aes(x = exon_classif, fill = intra_priming), 
                   stat = "count", color = "black", position = "dodge", width = 0.8) +
          labs(x = "Exon classification", y = "No. of transcripts") +
          theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
          scale_fill_manual("Intra-priming prediction", 
                            values = c("violetred3", "gray60")) +   
          facet_grid(~structural_category)
        
        # percentage
        ip_percent <- ggplot(ip_byexons) +
          geom_bar(aes(x = exon_classif, y = percent, 
                       fill = structural_category, alpha = intra_priming),
                   stat = "identity", position = "stack", 
                   width = 0.8, color = "black") +
          labs(x = "Exon classification", y = "% transcripts") +
          theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
          scale_fill_manual(guide = "none",
                            values = cat.palette) + 
          scale_alpha_manual("Intra-priming", values = c(0.3, 1)) +
          scale_y_continuous(labels = scales::percent_format()) +
          facet_grid(~structural_category)
        
#### END intra-priming plots #### 
    



####-------------------------- GENERATE PDF REPORT ------------------------####
        

