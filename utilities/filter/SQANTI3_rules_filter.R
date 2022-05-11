#!/bin/bash Rscript

##########################################################
#############       SQANTI3 RULES FILTER      ############
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
library(jsonlite)
library(magrittr)
library(dplyr)

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

### read files
message("-------------------------------------------------")
message("\n \t Reading classification file")
message("\n--------------------------------------------------")
classif <- read.table(classif_file, sep="\t", header = T, as.is = T)
message("-------------------------------------------------")
message("\n \t Reading JSON file with rules to filter")
message("\n--------------------------------------------------")
json_df <- jsonlite::fromJSON(txt = json_file ,simplifyDataFrame = T)

### json list will be transformed into a list of data frames with the rules. Rules DF have 4 columns:
### structural_category {any category + rest } ; 
### column { any SQ3 column name} ; 
### type {Category, Min_Threshold and Max_Threshold} ;
######### Ranges will be converted into 2 rules, a min_threshold and a max_threshold
######### "Category" rules will accept arrays in the json file and it will also deal with the TRUE or FALSE columns, treating them as characters
### rule {the value that will be used to filter}
###########################################################
### It's possible to define several independent rules for the same structural category, that's why I can't use structural categories to iterate the json list

rules_list <- list() 
list_idx <- 0
names_rules_list <- c()
for (sc in names(json_df)){
  for (j in c(1:dim(json_df[[sc]])[1])){
    rules <- json_df[[sc]][j,]
    names(rules) <- names(json_df[[sc]])
    rules_table <- data.frame()
    count=0
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
            rules_table[count,"rule"]=tolower(r[c])
            if (c != length(r)){
              count=count+1
            }
          }
        }
      }else if(is.numeric(r)){
        rules_table[count,"structural_category"]=sc
        rules_table[count,"column"]=col_name
        rules_table[count,"type"]="Min_Threshold"
        rules_table[count,"rule"]=r
      }else if(is.character(r)){
        rules_table[count,"structural_category"]=sc
        rules_table[count,"column"]=col_name
        rules_table[count,"type"]="Category"
        rules_table[count,"rule"]=tolower(r)
      }
    }
    list_idx=list_idx+1
    rules_table <- rules_table[!is.na(rules_table$rule),]
    names_rules_list <- c(names_rules_list, sc)
    rules_list[[list_idx]]<- rules_table
  }
}


names(rules_list) <- names_rules_list
## This is the actual function that will classify transcripts as Isoform or Artifacts

apply_rules <- function(isoform_info){
  sc = as.character(isoform_info["structural_category"])
  final_is_isoform=TRUE
  # detect if there is any specific rule for the SC of the isoform
  if (sc %in% names(json_df)) {
    final_is_isoform=FALSE
    # iterate all the independent rules for a certain SC
    for (p in which(names(json_df)==sc)){
      is_isoform=TRUE
      rules <- rules_list[[p]] 
      # iterate through the rules defined
      for (i in c(1:length(rules$rule))){
        if ( ! is.na(isoform_info[rules[i, "column"]])){ # if NA in the field, rule doesn't apply
          if (rules[i, "type"] == "Min_Threshold"){
            if (as.numeric(isoform_info[rules[i, "column"]]) < as.numeric(rules[i, "rule"])){
              is_isoform=FALSE
            }
          }else if (rules[i, "type"] == "Max_Threshold"){
            if (as.numeric(isoform_info[rules[i, "column"]]) > as.numeric(rules[i, "rule"])){
              is_isoform=FALSE
            }
          }else if (rules[i, "type"] == "Category"){
            cat_rules <- rules[rules$column == rules[i, "column"], ]
            if ( ! tolower(isoform_info[rules[i, "column"]]) %in% cat_rules[,"rule"]){
              is_isoform=FALSE
            }
          }
        }
      }
      final_is_isoform=final_is_isoform | is_isoform
    }
  # the isoform has a SC different, if will be evaluated with the rules of "rest" (if any was defined)  
  }else if ("rest" %in% names(json_df)){
    final_is_isoform=FALSE
    for (p in which(names(json_df)=="rest")){
      is_isoform=TRUE
      rules <- rules_list[[p]]
      for (i in c(1:length(rules$rule))){
        if (! is.na(isoform_info[rules[i, "column"]])){
          if (rules[i, "type"] == "Min_Threshold"){
            if (as.numeric(isoform_info[rules[i, "column"]]) < as.numeric(rules[i, "rule"])){
              is_isoform=FALSE
            }
          }else if (rules[i, "type"] == "Max_Threshold"){
            if (as.numeric(isoform_info[rules[i, "column"]]) > as.numeric(rules[i, "rule"])){
              is_isoform=FALSE
            }
          }else if (rules[i, "type"] == "Category"){
            cat_rules <- rules[rules$column == rules[i, "column"], ]
            if ( ! tolower(isoform_info[rules[i, "column"]]) %in% cat_rules[,"rule"]){
              is_isoform=FALSE
            }
          }
        }
      }
      final_is_isoform=final_is_isoform | is_isoform
    }
  }
  if (final_is_isoform){
   return("Isoform") 
  }else{
   return("Artifact")
  }
}

### This function will apply the same filtering process but it will return the reasons for filtering a transcript. Use only on artifacts
get_reasons <- function(isoform_info){
  sc = as.character(isoform_info["structural_category"])
  reasons <- c()
  # detect if there is any specific rule for the SC of the isoform
  if (sc %in% names(json_df)) {
    # iterate all the independent rules for a certain SC
    for (p in which(names(json_df)==sc)){
      rules <- rules_list[[p]] 
      # iterate through the rules defined
      for (i in c(1:length(rules$rule))){
        if ( ! is.na(isoform_info[rules[i, "column"]])){ # if NA in the field, rule doesn't apply
          if (rules[i, "type"] == "Min_Threshold"){
            if (as.numeric(isoform_info[rules[i, "column"]]) < as.numeric(rules[i, "rule"])){
              reason_line <- paste("Low", rules[i, "column"])
              reasons <- c(reasons, reason_line)
            }
          }else if (rules[i, "type"] == "Max_Threshold"){
            if (as.numeric(isoform_info[rules[i, "column"]]) > as.numeric(rules[i, "rule"])){
              reason_line <- paste("High", rules[i, "column"])
              reasons <- c(reasons, reason_line)            
              }
          }else if (rules[i, "type"] == "Category"){
            cat_rules <- rules[rules$column == rules[i, "column"], ]
            if ( ! tolower(isoform_info[rules[i, "column"]]) %in% cat_rules[,"rule"]){
              reason_line <- paste("Out", rules[i, "column"])
              reasons <- c(reasons, reason_line)            
              }
          }
        }
      }
    }
    # the isoform has a SC different, if will be evaluated with the rules of "rest" (if any was defined)  
  }else if ("rest" %in% names(json_df)){
    for (p in which(names(json_df)=="rest")){
      rules <- rules_list[[p]]
      for (i in c(1:length(rules$rule))){
        if (! is.na(isoform_info[rules[i, "column"]])){
          if (rules[i, "type"] == "Min_Threshold"){
            if (as.numeric(isoform_info[rules[i, "column"]]) < as.numeric(rules[i, "rule"])){
              reason_line <- paste("Low", rules[i, "column"])
              reasons <- c(reasons, reason_line)
            }
          }else if (rules[i, "type"] == "Max_Threshold"){
            if (as.numeric(isoform_info[rules[i, "column"]]) > as.numeric(rules[i, "rule"])){
              reason_line <- paste("High", rules[i, "column"])
              reasons <- c(reasons, reason_line)            
            }
          }else if (rules[i, "type"] == "Category"){
            cat_rules <- rules[rules$column == rules[i, "column"], ]
            if ( ! tolower(isoform_info[rules[i, "column"]]) %in% cat_rules[,"rule"]){
              reason_line <- paste("Out", rules[i, "column"])
              reasons <- c(reasons, reason_line)            
            }
          }
        }
      }
    }
  }
 num_reasons <- length(unique(reasons))
 final_df <- data.frame(isoform=rep(isoform_info["isoform"], num_reasons),
                        structural_category=rep(isoform_info["structural_category"],num_reasons), 
                        reasons=unique(reasons))
 return(final_df)
}



message("-------------------------------------------------")
message("\n \t Performing filtering")
message("\n--------------------------------------------------")


classif$filter_result <- apply(classif, 1, apply_rules)

inclusion_list <- classif[classif$filter_result == "Isoform", "isoform"]

artifacts_classif <- classif[classif$filter_result == "Artifact", ]
reasons_df <- apply(artifacts_classif,1, get_reasons) %>% bind_rows()

message("-------------------------------------------------")
message("\n \t Writting results")
message("\n--------------------------------------------------")

write.table(classif, file=paste0(opt$dir, "/", opt$output, "_RulesFilter_result_classification.txt"),
            quote = FALSE, col.names = TRUE, sep ='\t', row.names = FALSE)


write.table(inclusion_list, file = paste0(opt$dir, "/", opt$output, "_inclusion-list.txt"),
            quote = FALSE, col.names = FALSE, sep ='\t', row.names = FALSE)

write.table(reasons_df, file = paste0(opt$dir, "/", opt$output, "_filtering_reasons.txt"),
            quote = FALSE, col.names = TRUE, sep ='\t', row.names = FALSE)




