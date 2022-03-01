#!/bin/bash Rscript

##########################################################
#########     SQANTI3 MACHINE LEARNING FILTER    #########
##########################################################
#
# Authors: Francisco J. Pardo-Palacios
# Contact: f.pardo.palacios@gmail.com
# Affiliation: Institute for Integrative Systems Biology, CSIC, Valencia, Spain
#
# Last updated: Sept/27/2021
#
# Rules filter will take as input a JSON file where the user will set the rules for the filtering.
# It will be a JSON with nested objects so you can define any column you want to use for filtering.
# Filtering values will ALWAYS be in "positive", so the values set in the JSON will be the attributes that the filter will keep

library(optparse)
library(tidyverse)
library(jsonlite)

### Define arguments
option_list = list(
  optparse::make_option(c("-c","--sqanti_classif"), type="character", default = NULL, 
                        help="SQANTI3 classification output file." ),
  optparse::make_option(c("-j","--json_filter"), type="character",default = NULL,
                        help="Path to JSON file containing the rules for filtering isoforms. 
              If it is not provided, it will be used the default one stored in utilities folder."),
  optparse::make_option(c("-o","--output"), type="character", default = "SQANTI3", 
                        help="Output classification file prefix."),
  optparse::make_option(c("-d","--dir"), type="character", 
                        help="Output directory.")
)


### Parse and read arguments

opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser)
classif_file = opt$sqanti_classif
json_file = opt$json_filter

classif <- read.table(classif_file, sep="\t", header = T, as.is = T)

json_df <- jsonlite::fromJSON(txt = "test.json" ,simplifyDataFrame = T)

### read json and convert into data frame with 4 columns:
### structural_category {any category + rest } ; 
### column { any SQ3 column name} ; 
### type {Logical, Category, Min_Threshold and Max_Threshold}
### rule {the value that will be use to filter}

message("-------------------------------------------------")
message("\n \t Reading JSON file with rules to filter")
message("\n--------------------------------------------------")


rules_table <- data.frame()
count=0
for (sc in names(json_df)){
  rules <- json_df[[as.character(sc)]]
  for (r in names(rules)){
    count=count+1
    col_name <- r
    r <- unlist(rules[[r]])
    if (length(r)>1){
      if (is.numeric(r)){
        lower_limit <- min(r)
        upper_limit <- max(r)
        rules_table[count,"structural_category"]=sc
        rules_table[count,"column"]=col_name
        rules_table[count,"type"]="Min_Threshold"
        rules_table[count,"rule"]=lower_limit
        count=count+1
        rules_table[count,"structural_category"]=sc
        rules_table[count,"column"]=col_name
        rules_table[count,"type"]= "Max_Threshold"
        rules_table[count,"rule"]= upper_limit
      }else{
        filter_type <- "Category"
        for (c in c(1:length(r))){
          rules_table[count,"structural_category"]=sc
          rules_table[count,"column"]=col_name
          rules_table[count,"type"]=filter_type
          rules_table[count,"rule"]=r[c]
          if (c != length(r)){
            count=count+1
          }
        }
      }
    }else if(r=="TRUE"|r=="FALSE"|r=="True"|r=="False"){
      rules_table[count,"structural_category"]=sc
      rules_table[count,"column"]=col_name
      rules_table[count,"type"]="Logical"
      rules_table[count,"rule"]=as.logical(r)
    }else if(is.numeric(r)){
      rules_table[count,"structural_category"]=sc
      rules_table[count,"column"]=col_name
      rules_table[count,"type"]="Min_Threshold"
      rules_table[count,"rule"]=r
    }else if(is.character(r)){
      rules_table[count,"structural_category"]=sc
      rules_table[count,"column"]=col_name
      rules_table[count,"type"]="Category"
      rules_table[count,"rule"]=r
    }
  }
}

## Define the actual function that will classify transcripts as Isoform or Artifacts

apply_rules <- function(isoform_info){
  o="Isoform"
  sc = as.character(isoform_info["structural_category"])
  if (sc %in% rules_table[,"structural_category"]){
      rules <- rules_table[rules_table$structural_category == sc , ]
      for (i in c(1:length(rules$rule))){
        if ( ! is.na(isoform_info[rules[i, "column"]])){
           if (rules[i, "type"] == "Min_Threshold"){
                if (as.numeric(isoform_info[rules[i, "column"]]) < as.numeric(rules[i, "rule"])){
                    o="Artifact"
                }
           }else if (rules[i, "type"] == "Max_Threshold"){
                 if (as.numeric(isoform_info[rules[i, "column"]]) > as.numeric(rules[i, "rule"])){
                    o="Artifact"
                }
           }else if (rules[i, "type"] == "Logical"){
                 if (as.logical(isoform_info[rules[i, "column"]]) != rules[i, "rule"]){
                    o="Artifact"
                 }
           }else if (rules[i, "type"] == "Category"){
                 cat_rules <- rules[rules$column == rules[i, "column"], ]
                 if ( ! as.character(isoform_info[rules[i, "column"]]) %in% cat_rules[,"rule"]){
                    o="Artifact"
                 }
           }
        }
      }    
  } else {
    rules <- rules_table[rules_table$structural_category == "rest", ]
    for (i in c(1:length(rules$rule))){
      if (! is.na(isoform_info[rules[i, "column"]])){
           if (rules[i, "type"] == "Min_Threshold"){
                if (as.numeric(isoform_info[rules[i, "column"]]) < as.numeric(rules[i, "rule"])){
                    o="Artifact"
                }
           }else if (rules[i, "type"] == "Max_Threshold"){
                 if (as.numeric(isoform_info[rules[i, "column"]]) > as.numeric(rules[i, "rule"])){
                    o="Artifact"
                }
           }else if (rules[i, "type"] == "Logical"){
                 if (as.logical(isoform_info[rules[i, "column"]]) != rules[i, "rule"]){
                    o="Artifact"
                 }
           }else if (rules[i, "type"] == "Category"){
                 cat_rules <- rules[rules$column == rules[i, "column"], ]
                 if ( ! as.character(isoform_info[rules[i, "column"]]) %in% cat_rules[,"rule"]){
                    o="Artifact"
                 }
           }
      }
    }
  }
  return(o)
}

message("-------------------------------------------------")
message("\n \t Performing filtering")
message("\n--------------------------------------------------")


classif$Rules_filter <- apply(classif, 1, apply_rules)

inclusion_list <- classif[classif$Rules_filter == "Isoform", "isoform"]

message("-------------------------------------------------")
message("\n \t Writting results")
message("\n--------------------------------------------------")

write.table(classif, file=paste0(opt$dir, "/", opt$output, "_RulesFilter_result_classification.txt"),
            quote = FALSE, col.names = TRUE, sep ='\t', row.names = FALSE)


write.table(inclusion_list, file = paste0(opt$dir, "/", opt$output, "_inclusion-list.txt"),
            quote = FALSE, col.names = FALSE, sep ='\t', row.names = FALSE)



