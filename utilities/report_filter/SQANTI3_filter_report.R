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
# Last updated: 23/August/2023
#
#-------------------------------------------------------------------------------



####------------------------- FILTER REPORT ARG LIST -----------------------####

#### Define script arguments ####

option_list <- list(
  optparse::make_option(c("-d","--dir"), type = "character", 
                        help = "Output/input directory - must be same as SQ3 filter."),
  optparse::make_option(c("-o","--output"), type = "character", default = "SQANTI3", 
                        help = "Output report file prefix."),
  optparse::make_option(c("-u", "--utilities_path"), type = "character",
                        help = "Full path to SQANTI3/utilities folder."),
  optparse::make_option(c("-f", "--filter_type"), type = "character", default = "ml",
                        help = "Type of SQ3 filter that was run to generate input.
                        Must be one of the following: 'ml' or 'rules'.
                        Default: 'ml'.")
)

# Parse arguments
opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser)




####------------------------ INPUTS & PLOT THEME ---------------------------####

#### Initialize script and load input data ####
if (opt$filter_type == "ml"){
  message("\n-------------------------------------------------")
  message("\n \t SQANTI3 Machine Learning filter report")
  message("\n--------------------------------------------------")
}else if (opt$filter_type == "rules"){
  message("\n-------------------------------------------------")
  message("\n \t SQANTI3 Rules filter report")
  message("\n--------------------------------------------------")
}


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
  
  # Detect path and load parameters
  message("\nReading ML filter parameters...")
  
  which_params <- stringr::str_detect(paths, "params")
  path_params <- paths[which_params]
  params <- readr::read_tsv(path_params, col_names = c("parameter", "value"))
  
  # Detect path and load model performance stats files
  message("\nReading ML performance statistics...")
  
    # statistics
    which_stats <- stringr::str_detect(paths, "testSet_stats")
    path_stats <- paths[which_stats]
    stats <- readr::read_tsv(path_stats, col_names = c("metric", "value"))
    
    # confusion matrix
    which_conf <- stringr::str_detect(paths, "confusionMatrix")
    path_conf <- paths[which_conf]
    conf <- readr::read_tsv(path_conf)
    
  
}else if (opt$filter_type == "rules"){
  # Detect path and load ML output classification
  message("\nReading Rules result classification table...")
  
  which_classif <- stringr::str_detect(paths, "RulesFilter_result_classification")
  path_classif <- paths[which_classif]
  classif <- readr::read_tsv(path_classif, 
                             col_types = readr::cols(exons = readr::col_integer(),
                                                     ref_exons = readr::col_integer()))
  
}



# Format category factor for classification file
classif <- classif %>% 
  dplyr::mutate(structural_category = factor(structural_category) %>% 
                  forcats::fct_recode(FSM = "full-splice_match",
                                      ISM = "incomplete-splice_match",
                                      NIC = "novel_in_catalog",
                                      NNC = "novel_not_in_catalog",
                                      Intergenic = "intergenic",
                                      Antisense = "antisense",
                                      Genic = "genic",
                                      Fusion = "fusion",
                                      `Genic \nintron` = "genic_intron") %>% 
                  forcats::fct_relevel(c("FSM", "ISM", "NIC", "NNC", 
                                         "Genic", "Antisense", "Fusion", 
                                         "Intergenic", "Genic \nintron"))) 


#### Define plot theme ####

# Load ggplot2
require(ggplot2)

# Install RColorConesa if not available
pkg <- installed.packages() %>% rownames

if(!("RColorConesa" %in% pkg)){
  suppressMessages(install.packages("RColorConesa"))
}

# Set theme parameters (from SQANTI3_report.R)
sq_theme <- theme_classic(base_family = "Helvetica") +
  theme(plot.title = element_text(lineheight=.4, size=15, hjust = 0.5)) +
  theme(plot.margin = unit(c(2.5,1,1,1), "cm")) +
  theme(axis.line.x = element_line(color="black", linewidth = 0.4),
        axis.line.y = element_line(color="black", linewidth = 0.4)) +
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
    
message("\nGenerating common filter plots...\n")

# Gene-level summary tables
gene_artifacts <- classif %>% 
  dplyr::select(associated_gene, filter_result) %>% 
  dplyr::mutate(gene_type = ifelse(stringr::str_detect(associated_gene, "novel"),
                                   yes = "Novel", no = "Annotated")) %>% 
  dplyr::group_by(associated_gene) %>% 
  dplyr::reframe(all_artifacts = all(filter_result == "Artifact"),
                   gene_type = gene_type,
                   filter_result = filter_result)

bygene_summary <- gene_artifacts %>% 
  dplyr::select(!filter_result) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(gene_type) %>% 
  dplyr::add_tally(name = "total_genes") %>% 
  dplyr::group_by(gene_type, all_artifacts) %>% 
  dplyr::reframe(isoform_no = dplyr::n(),
                 total_genes = total_genes) %>%  
  dplyr::distinct()

# long-format table for grid.table()
bygene_summary_long <- bygene_summary %>% 
  tidyr::pivot_wider(names_from = "all_artifacts", 
                     names_prefix = "all_artifacts_", 
                     values_from = "isoform_no") %>% 
  dplyr::select(gene_type, total_genes, all_artifacts_TRUE)



# Isoforms and artifacts by category
category_summary <- classif %>% 
  dplyr::group_by(structural_category, filter_result) %>% 
  dplyr::summarize(n = dplyr::n()) %>%
  dplyr::mutate(percent = n/sum(n))

# long-format table for grid.table()
category_summary_long <- category_summary %>% 
  dplyr::select(-percent) %>% 
  tidyr::pivot_wider(names_from = filter_result, values_from = n)



      ### Create table objects for report ###
      # Total genes and artifacts 
      gene_count <- classif$associated_gene %>% unique %>% length
      tr_count <- classif %>% 
        dplyr::group_by(filter_result) %>% 
        dplyr::summarize(n = dplyr::n()) %>% 
        tibble::deframe()
      
      iso_count <- tr_count["Isoform"]
      iso_pcnt <- iso_count*100/sum(tr_count)
      artif_count <- tr_count["Artifact"]
      artif_pcnt <- artif_count*100/sum(tr_count) %>% round
      
      
      sentence <- paste0("Total Genes: ", gene_count, "\n\n", 
                         "Total Transcripts: ", sum(tr_count), "\n",
                         "- Isoforms: ", 
                         iso_count, " (", round(iso_pcnt), "%)", "\n",
                         "- Artifacts: ", 
                         artif_count, " (", round(artif_pcnt), "%)")
      
      summary_title <- grid::textGrob(sentence, 
                                      gp = grid::gpar(fontface = "italic", 
                                                      fontsize = 17), vjust = 0)
      
      # Total artifacts/isoforms by category
      gene_table <- gridExtra::tableGrob(bygene_summary_long, rows = NULL,
                                         cols = c("Gene category", "Gene no.", 
                                                  "No. of genes with \nartifacts only"))
      
      # Total genes and genes with all artifacts
      cat_table <- gridExtra::tableGrob(category_summary_long, rows = NULL,
                                        cols = c("Structural category",
                                                 "Artifact no.", "Isoform no."))


    
      ### Create summary plots for report ###
      # Plot totals
      cat_totals <- ggplot(category_summary, 
                            aes(x = structural_category, y = n)) + 
        ggtitle("Total isoforms and artifacts by category") +
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
        ggtitle("% isoforms and artifacts by category") +
        geom_bar(aes(alpha = filter_result, fill = structural_category), 
                  stat = "identity",
                  width = 0.8, color = "black", position = "stack") +
        labs(x = "Structural category", y = "Transcript %") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        scale_y_continuous(labels = scales::percent_format()) +
        scale_fill_manual("", values = cat.palette, guide = "none") +
        scale_alpha_manual("Filter result", values = c(0.3, 1))

      
      
# Number of isoforms per gene before and after filtering
# calculate before
isoforms_before <- classif %>% 
  dplyr::select(isoform, associated_gene) %>% 
  dplyr::group_by(associated_gene) %>% 
  dplyr::summarize(isoform_no = dplyr::n()) %>% 
  dplyr::mutate(isoform_factor = cut(isoform_no, c(0, 1, 3, 5, max(.$isoform_no)),
                                     labels = c("1", "2-3", "4-5", ">=6")))

# calculate after
isoforms_after <- classif %>% 
  dplyr::filter(filter_result == "Isoform") %>%  
  dplyr::group_by(associated_gene) %>% 
  dplyr::summarize(isoform_no = dplyr::n()) %>% 
  dplyr::mutate(isoform_factor = cut(isoform_no, c(0, 1, 3, 5, max(.$isoform_no)),
                                     labels = c("1", "2-3", "4-5", ">=6")))

# join both
isoforms_per_gene <- dplyr::bind_rows(list(Before = isoforms_before, 
                                           After = isoforms_after),
                                      .id = "filter") %>% 
  dplyr::mutate(filter = factor(filter) %>% 
                  forcats::fct_relevel(c("Before", "After")))

# compute percentage
isoforms_per_gene <- isoforms_per_gene %>% 
  dplyr::group_by(filter, isoform_factor) %>% 
  dplyr::summarize(gene_no = dplyr::n()) %>% 
  dplyr::mutate(percent = gene_no/sum(gene_no))

      # plot
      iso_p_gene <- ggplot(isoforms_per_gene) +
        ggtitle("Number of isoforms per gene") +
        geom_bar(aes(x = isoform_factor, y = percent, fill = filter),
                 stat = "identity", position = "dodge",
                 color = "black", width = 0.8) +
        labs(x = "Isoforms per gene", y = "% Genes") +
        scale_fill_manual("Filter", values = c("tomato3", "#FFD4BF")) +
        scale_y_continuous(labels = scales::percent_format())

      
      
# Redundancy before and after filtering
# fsm redundancy
fsm_before <- classif %>% 
  dplyr::filter(structural_category == "FSM")

fsm_after <- classif %>% 
  dplyr::filter(structural_category == "FSM" &
                  filter_result == "Isoform")

fsm <- dplyr::bind_rows(list(Before = fsm_before, After = fsm_after),
                        .id = "filter") %>% 
  dplyr::mutate(filter = factor(filter) %>% 
                  forcats::fct_relevel(c("Before", "After"))) %>% 
  dplyr::select(isoform, associated_gene, associated_transcript, 
                structural_category, filter)

fsm_redund.df <- fsm %>% 
  dplyr::group_by(filter, associated_transcript) %>% 
  dplyr::summarize(fsm_no = dplyr::n()) %>% 
  dplyr::filter(fsm_no > 0) %>% 
  dplyr::mutate(fsm_type = as.factor(ifelse(fsm_no == 1, 
                                            yes = "Unique", no = "Multiple")) %>% 
                  forcats::fct_relevel(c("Unique", "Multiple")),
                fsm_bin = ifelse(fsm_no <= 4, 
                                 yes = as.character(fsm_no), no = ">4") %>% 
                  ordered(levels = c("1", "2", "3", "4", ">4")))
      
      # plot general reference transcript complexity
      ref_complexity <- ggplot(fsm_redund.df) +
        ggtitle("Reference transcript complexity", 
                subtitle = "No. of reference transcripts represented by FSM") +
        geom_bar(aes(x = filter, fill = fsm_type), 
                 stat = "count", position = "stack",
                 color = "black", width = 0.7) +
        geom_text(aes(x = filter,
                      label = after_stat(count)), 
                  stat = "count", vjust = -1) +
        scale_y_continuous(breaks = scales::pretty_breaks(6)) +
        RColorConesa::scale_fill_conesa("FSM per \nreference ID", palette = "nature",
                                        continuous = FALSE, reverse = FALSE) +
        theme(plot.subtitle = element_text(hjust = 0.5, size = 12,
                                           face = "italic")) +
        xlab("Filter") +
        ylab("Reference transcfript no.") 
        
      
      # plot FSM redundancy
      fsm_redund <- ggplot(fsm_redund.df) + 
        ggtitle("FSM redundancy") +
        geom_bar(aes(x = fsm_type, fill = fsm_bin), color = "black", width = 0.5) +
        geom_text(aes(x = fsm_type, label = after_stat(count)), stat = "count", vjust = -1) +
        scale_y_continuous(breaks = scales::pretty_breaks(6)) +
        scale_fill_brewer("Total FSM \nper reference ID", palette = "Blues") +
        xlab("FSM per reference transcript") +
        ylab("Reference transcfript no.") +
        facet_grid(~filter)


# ism redundancy
ism_before <- classif %>% 
  dplyr::filter(structural_category == "ISM")

ism_after <- classif %>% 
  dplyr::filter(structural_category == "ISM" &
                  filter_result == "Isoform")

ism <- dplyr::bind_rows(list(Before = ism_before, After = ism_after),
                        .id = "filter") %>% 
  dplyr::mutate(filter = factor(filter) %>% 
                  forcats::fct_relevel(c("Before", "After"))) %>% 
  dplyr::select(isoform, associated_gene, associated_transcript, 
                structural_category, filter)

ism_redund.df <- ism %>% 
  dplyr::group_by(filter, associated_transcript) %>% 
  dplyr::summarize(total_ism = dplyr::n()) %>% 
  dplyr::filter(total_ism > 0) %>% 
  dplyr::mutate(ism_type = ifelse(total_ism == 1, 
                                  yes = "Unique", no = "Multiple") %>% 
                  forcats::fct_relevel(c("Unique", "Multiple")),
                ism_bin = ifelse(total_ism <= 4, 
                                 yes = as.character(total_ism), no = ">4") %>% 
                  ordered(levels = c("1", "2", "3", "4", ">4")))

      # plot FSM redundancy
      ism_redund <- ggplot(ism_redund.df) + 
        ggtitle("ISM redundancy") +
        geom_bar(aes(x = ism_type, fill = ism_bin), color = "black", width = 0.5) +
        geom_text(aes(x = ism_type, label = after_stat(count)), stat = "count", vjust = -1) +
        scale_y_continuous(breaks = scales::pretty_breaks(6)) +
        scale_fill_brewer("Total ISM \nper reference ID", palette = "Oranges") +
        xlab("ISM per reference transcript") +
        ylab("Reference transcfripts") +
        facet_grid(~filter)

      
# combined redundancy (FSM + ISM)
comb_fsm_ism <- dplyr::bind_rows(fsm, ism)

comb_redund.df <- comb_fsm_ism %>% 
  dplyr::group_by(filter, associated_transcript) %>% 
  dplyr::summarize(total_tr = dplyr::n()) %>% 
  dplyr::filter(total_tr > 0) %>% 
  dplyr::mutate(tr_type = ifelse(total_tr == 1, 
                                  yes = "Unique", no = "Multiple") %>% 
                  forcats::fct_relevel(c("Unique", "Multiple")),
                tr_bin = ifelse(total_tr <= 4, 
                                 yes = as.character(total_tr), no = ">4") %>% 
                  ordered(levels = c("1", "2", "3", "4", ">4")))
      
        # plot combined redundancy
        comb_redund <- ggplot(comb_redund.df) + 
          ggtitle("FSM+ISM redundancy") +
          geom_bar(aes(x = tr_type, fill = tr_bin), color = "black", width = 0.5) +
          geom_text(aes(x = tr_type, label = after_stat(count)), stat = "count", vjust = -1) +
          scale_y_continuous(breaks = scales::pretty_breaks(6)) +
          scale_fill_brewer("Total FSM+ISM \nper reference ID", palette = "Greens") +
          xlab("FSM+ISM per reference transcript") +
          ylab("Reference transcfript no.") +
          facet_grid(~filter)
        
        

#### END common filter plots ####
      
      
      
#### ML filter plots ####
if(opt$filter_type == "ml"){
    
    message("\nGenerating machine learning filter plots...\n")
    
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
          labs(x = "Structural category", y = "Transcript no.") +
          theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
          scale_fill_manual("Reason", values = c("lightcoral", "lightgoldenrod1",
                                                 "aquamarine3"))
        
        # percentage
        artifact_percent <- ggplot(artifact_summary) +
          ggtitle("Reason to flag transcripts as artifacts, by category") +
          geom_bar(aes(x = structural_category, y = percent, fill = reason), 
                   stat = "identity", position = "stack", 
                   color = "black", width = 0.8) +
          labs(x = "Structural category", y = "Transcript %") +
          theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
          scale_fill_manual("Reason", values = c("lightcoral", "lightgoldenrod1",
                                                 "aquamarine3")) +
          scale_y_continuous(labels = scales::percent_format())
        
    
    ## Model performance: summary table of test set stats metrics
    performance_title <- grid::textGrob("Classification model performance on test set", 
                                  gp = grid::gpar(fontface = "bold", 
                                                        fontsize = 15), vjust = 0)
    
      # stats with changed number formatting
      # requires handling p-value separately (values too low)
      accpvalue <- stats %>% 
        dplyr::filter(metric == "AccuracyPValue")
      pvalue_table <- gridExtra::tableGrob(accpvalue, rows = NULL,
                                           cols = NULL)
      
      stats <- stats %>% 
        dplyr::mutate(value = tibble::num(stats$value, sigfig = 3, 
                                          notation = "dec")) %>% 
        dplyr::filter(metric != "AccuracyPValue")
      stats_title <- grid::textGrob("Performance metrics", 
                                    gp = grid::gpar(fontsize = 14), vjust = 0)
      stats_table <- gridExtra::tableGrob(stats, rows = NULL,
                                          cols = c("Metric", "Value"))
      
    
      # confusion matrix
      conf_title <- grid::textGrob("Confusion matrix",
                                   gp = grid::gpar(fontsize = 14), vjust = 0)
      conf_table <- gridExtra::tableGrob(conf, rows = NULL)
      gconf <- gridExtra::arrangeGrob(conf_title, conf_table,
                                       layout_matrix = cbind(c(1,2,3)),
                                       heights = c(0.2,0.5,1))
    
            
    ## Variable importance in classifier
    var_imp <- ggplot(imp) +
      ggtitle("Variable importance in Random Forest classifier") +
      geom_bar(aes(x = variable, y = importance), stat = "identity",
               width = 0.7, color = "black", fill = "cadetblue3") +
      labs(x = "SQANTI3 variables", y = "Importance") +
      coord_flip()
    
    
    ## Probability distribution obtained from model
    probabilities <- classif %>% 
      dplyr::select(isoform, POS_MLprob) %>% 
      dplyr::filter(!is.na(POS_MLprob))
    
    prob_dens <- ggplot(probabilities) +
      ggtitle(label = "",
        subtitle = paste("Transcripts classified:", nrow(probabilities))) +
      geom_density(aes(x = POS_MLprob)) +
      labs(x = "Positive/Isoform probability obtained from ML filter classifier (POS_MLprob)",
           y = "Density")
    
    ## Variables used in ML: values for isoforms and artifacts by category
    source(paste0(opt$utilities_path, "/report_filter/compare_MLvariables.R"))
    
    var_compare <- purrr::map2(imp$variable, imp$importance,
                               ~compare_MLvariables(classif, .x, .y))
    
}
#### END ML filter plots ####

      
      
#### ML Intra-priming plots ####
if(opt$filter_type == "ml"){
  
  message("\nGenerating ML filter intra-priming plots...\n")
  
  # A% per category
  a_value <- dplyr::filter(params, parameter == "Intrapriming") %>% 
    dplyr::select(value) %>% tibble::deframe() %>% 
    as.numeric()
  
  # plot
  a_percent <- ggplot(classif) +
    ggtitle("A % by category", 
            subtitle = "Red line indicates threshold employed in ML filter") +
    geom_boxplot(aes(x = structural_category, y = perc_A_downstream_TTS),
                 width = 0.5, fill = "navajowhite1") +
    geom_hline(aes(yintercept = a_value), 
               color = "firebrick3", size = 1, linetype = "dashed") +
    xlab("Structural category") +
    ylab("% A's downstream of TTS")
  
  
  
  # Intra-priming by category and exon number
  ip_summary <- classif %>% 
    dplyr::select(structural_category, exons, intra_priming) %>% 
    dplyr::mutate(exon_classif = dplyr::if_else(exons > 1, 
                                                true = "Multi-exon", 
                                                false = "Mono-exon"))
  
  ip_general <- ip_summary %>% 
    dplyr::group_by(structural_category, intra_priming) %>%
    dplyr::summarize(n = dplyr::n()) %>% 
    dplyr::mutate(percent = n/sum(n))
  
  ip_byexons <- ip_summary %>% 
    dplyr::group_by(structural_category, exon_classif, intra_priming) %>% 
    dplyr::summarize(n = dplyr::n()) %>% 
    dplyr::mutate(percent = n/sum(n))
  
  
  # Intra-priming by structural_category
  # totals
  ip_totals <- ggplot(ip_summary) +
    ggtitle("Isoforms flagged as intra-priming, by category") +
    geom_bar(aes(x = structural_category, fill = intra_priming), 
             stat = "count", color = "black", position = "dodge", width = 0.8) +
    labs(x = "Exon classification", y = "Transcript no.") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    scale_fill_manual("Intra-priming prediction", 
                      values = c("#74CDF0", "#9F7BB8")) +
    scale_y_continuous(breaks = scales::pretty_breaks(6))
  
  # percentage
  ip_percent <- ggplot(ip_general) +
    ggtitle("Isoforms flagged as intra-priming, by category (%)") +
    geom_bar(aes(x = structural_category, y = percent, 
                 fill = structural_category, alpha = intra_priming),
             stat = "identity", position = "stack", 
             width = 0.8, color = "black") +
    labs(x = "Exon classification", y = "Transcript %") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    scale_fill_manual(guide = "none",
                      values = cat.palette) + 
    scale_alpha_manual("Intra-priming", values = c(0.3, 1)) +
    scale_y_continuous(labels = scales::percent_format())
  
  
  # Intra-priming by structural_category and exon_classif
  # totals
  ipex_totals <- ggplot(ip_summary) +
    ggtitle("Isoforms flagged as intra-priming, by exon number") +
    geom_bar(aes(x = exon_classif, fill = intra_priming), 
             stat = "count", color = "black", position = "dodge", width = 0.8) +
    labs(x = "Exon classification", y = "Transcript no.") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    scale_fill_manual("Intra-priming prediction", 
                      values = c("#74CDF0", "#9F7BB8")) +   
    facet_grid(~structural_category)
  
  # percentage
  ipex_percent <- ggplot(ip_byexons) +
    ggtitle("Isoforms flagged as intra-priming, by exon number (%)") +
    geom_bar(aes(x = exon_classif, y = percent, 
                 fill = structural_category, alpha = intra_priming),
             stat = "identity", position = "stack", 
             width = 0.8, color = "black") +
    labs(x = "Exon classification", y = "Transcript %") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    scale_fill_manual(guide = "none",
                      values = cat.palette) + 
    scale_alpha_manual("Intra-priming", values = c(0.3, 1)) +
    scale_y_continuous(labels = scales::percent_format()) +
    facet_grid(~structural_category)
  
}

#### END intra-priming plots #### 
        
####-------------------------- GENERATE PDF REPORT ------------------------####

message("\nWriting report plots to PDF file...")
        
# Create report file name
pdf_file <- paste0(opt$dir, "/", opt$output, "_SQANTI3_filter_report.pdf")

# Open file
pdf(file = pdf_file, width = 8, height = 7.5)

    # Print title and header
    grid::grid.newpage()
    cover <- grid::textGrob("SQANTI3 filter report",
                            gp = grid::gpar(fontface = "italic", 
                                            fontsize = 40, col = "orangered"))
    grid::grid.draw(cover)

    
# Create grid of tables
    
    gridExtra::grid.arrange(summary_title, gene_table, cat_table,
                            layout_matrix = cbind(c(1,2),c(1,4)))

# Print plots

    # Common filter plots
    print(cat_totals)
    print(cat_percent)
    print(iso_p_gene)
    print(ref_complexity)
    print(fsm_redund)
    print(ism_redund)
    print(comb_redund)
    
    # ML filter plots
    if(opt$filter_type == "ml"){
      
      # ML filter plots cover
      grid::grid.newpage()
      mlcover <- grid::textGrob("ML classifier performance report",
                              gp = grid::gpar(fontface = "italic", 
                                              fontsize = 30, col = "steelblue"))
      grid::grid.draw(mlcover)
      
      # grid of ML tables
      gridExtra::grid.arrange(performance_title, stats_title, 
                              stats_table, gconf, pvalue_table,
                              layout_matrix = cbind(c(1,2,3,5),c(1,4,4,4)),
                              heights = c(0.2,0.1,1,0.1))
      
      # plots
      print(prob_dens)
      print(var_imp)
      print(artifact_totals)
      print(artifact_percent)
      suppressWarnings(purrr::walk(var_compare, print))
    
      # Intra-priming plots (only available when ML activated)
        
      # ML filter plots cover
        grid::grid.newpage()
        ipcover <- grid::textGrob("Intra-primming filter report",
                                gp = grid::gpar(fontface = "italic", 
                                                fontsize = 30, col = "steelblue"))
        grid::grid.draw(ipcover)
      
      print(a_percent)
      print(ip_totals)
      print(ip_percent)
      print(ipex_totals)
      print(ipex_percent)
}
# Close file
dev.off()

        
