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
                        help = "Type of SQ3 filter that was run to generate input.")
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
  classif <- readr::read_tsv(path_classif)
  
  
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
                                      Genic_intron = "genic_intron"))



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
                Fusion="goldenrod1", Intergenic = "darksalmon", Genic_intron="#41B6C4")



#### Common filter plots ####
    


#### ML filter plots ####
if(opt$filter_type == "ml"){
  
    ## Variable importance in classifier
    var_imp <- ggplot(imp) +
      ggtitle("Variable importance in Random Forest classifier") +
      geom_bar(aes(x = variable, y = importance), stat = "identity",
               width = 0.7, color = "black", fill = "cadetblue3") +
      labs(x = "SQANTI3 variables", y = "Importance") +
      coord_flip()

    
}
#### END ML filter plots ####


##### Summary of results #####
# plot.compare is a function to evaluate the values of the features used in the 
# ML filter in Artifacts and Isoforms per structural category


plot.compare <- function (classification, myfeature, imp) {
  
  require(ggplot2)
  
  mycategories = c("full-splice_match", "incomplete-splice_match", "novel_in_catalog",
                   "novel_not_in_catalog", "intergenic", "fusion", "genic", 
                   "antisense", "genic_intron")
  mylabels = c("FSM", "ISM", "NIC", "NNC", "Inter", "Fus", "Gen", "Anti", "Intron")
  col <- which(colnames(classification) == myfeature)
  
  
  if (class(classification[,col]) != "factor") {
    if (myfeature != "count")  {
      dat <- data.frame(category = classification$structural_category, 
                        filter = classification$MLfilter_result, 
                        feature = classification[,col])
      dat2 <- with(dat, tapply(feature,list(category,filter), median))
    } else {
      dat <- data.frame(category = classification$structural_category, 
                        filter = classification$MLfilter_result)
      dat2 <- table(dat)
      dat2 <- dat2[mycategories,]
    }
    
    dat3 <- reshape::melt(dat2) ; names(dat3) <- c("category", "filter", "feature")
    dat3$category <- factor(dat3$category, levels = mycategories)
    dat3$filter <- factor(dat3$filter, levels = c("Isoform", "Artifact"))
    
    p <- ggplot(dat3, aes(x = category, y =  feature)) +
      geom_bar(
        aes(color = filter, fill = filter),
        stat = "identity", position = position_dodge(0.8)) +
      scale_color_manual(values = c("#0073C2FF", "#EFC000FF")) +
      scale_fill_manual(values = c("#0073C2FF", "#EFC000FF")) +
      ggtitle(paste(myfeature, "importance", imp, sep = " ")) + 
      scale_x_discrete(breaks = mycategories,
                       labels = mylabels)
    
  } else {   # code for categorical features
    col <- which(colnames(classification) == myfeature)
    dat <- data.frame(category = classification$structural_category, 
                      filter = classification$MLfilter_result, feature = classification[,col])
    pasted <- as.data.frame(table(apply(dat,1, paste, collapse = "#")))
    splitted <- strsplit(as.character(pasted$Var1), split = "#")
    dat2 <- data.frame(matrix(unlist(splitted), nrow = length(splitted), byrow = TRUE))
    names(dat2) <- c("category", "filter", "feature")
    dat3 <- data.frame(dat2, value = pasted$Freq)
    dat3$category <- factor(dat3$category, levels = mycategories)
    dat3$filter <- factor(dat3$filter, levels = c("Isoform", "Artifact"))
    dat3$labels <- mylabels [match(dat3$category,  mycategories)]
    dat3$labels <- factor(dat3$labels, levels = mylabels)
    
    p <- ggplot(dat3, aes(x = filter, y = value, fill = feature)) + 
      geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ labels) +
      ggtitle(paste(myfeature, "importance", imp, sep = " ")) 
  }
  
  return(p)
}

features.eval <- c(rownames(imp), "intra_priming")
importance <- c(round(imp[,1],2), "NA")

pdf(file = paste0(opt$output_directory, "/SQANTI3_ML_report.pdf"))
for (i in 1:length(features.eval)) {
  message(paste(i, ":", features.eval[i]))
  eval_p <- plot.compare(classification = d1, myfeature = features.eval[i], imp = importance[i])
  print(eval_p)
}
dev.off()


plot(density(log(d1$ratio_TSS[d1$structural_category == "incomplete-splice_match"])), col = "red")
lines(density(log(d1$ratio_TSS[d1$structural_category == "full-splice_match"])))
lines(density(log(d1$ratio_TSS[(d1$structural_category == "incomplete-splice_match" & d1$MLfilter_result == "Isoform")])), col = "blue")





