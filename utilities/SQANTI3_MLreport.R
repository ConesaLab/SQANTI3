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




plot(density(log(d1$ratio_TSS[d1$structural_category == "incomplete-splice_match"])), col = "red")
lines(density(log(d1$ratio_TSS[d1$structural_category == "full-splice_match"])))
lines(density(log(d1$ratio_TSS[(d1$structural_category == "incomplete-splice_match" & d1$MLfilter_result == "Isoform")])), col = "blue")





