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
# It will be a JSON with nested objects so you can define any column you want to use for filtering 

### Define arguments
option_list = list(
  optparse::make_option(c("-c","--sqanti_classif"), type="character", default = NULL, 
                        help="SQANTI classification output file." ),
  optparse::make_option(c("-j","--json_filter"), type="character",default = NULL,
                        help="Path to JSON file containing the rules for filtering isoforms. 
              If it is not provided, it will be used the default one stored in utilities folder.")
)



# testing
setwd("/home/fjpardo/Florida/hipergator3/github/SQANTI3/utilities/")
library(jsonlite)
library(tidyverse)
json_df <- jsonlite::fromJSON(txt = "test.json" ,simplifyDataFrame = T)

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
        filter_type <- "Max_Threshold"
        filter_rule <- upper_limit
      }else{
        filter_type <- "Category"
        for (c in length(r)-1){
          rules_table[count,"structural_category"]=sc
          rules_table[count,"column"]=col_name
          rules_table[count,"type"]=r[c]
          rules_table[count,"rule"]=r[c]
          count=count+1
        }
        filter_rule <- r[c+1]
      }
    }else if(r=="TRUE"|r=="FALSE"|r=="True"|r=="False"){
      filter_type <- "Logical"
      filter_rule <- as.logical(r)
    }else if(is.numeric(r)){
      filter_type <- "Min_Threshold"
      filter_rule <- r
    }else if(is.character(r)){
      filter_type <- "Category"
      filter_rule <- r
    }
    rules_table[count,"structural_category"]=sc
    rules_table[count,"column"]=col_name
    rules_table[count,"type"]=filter_type
    rules_table[count,"rule"]=filter_rule
  }
}

get_rules <- function(isoform_info){
  sc = isoform_info["structural_category"]
  if (sc %in% rules_table["structural_category"])){
      
  }else{
    sc="rest"
  }
}

filter_by_rules <- function(classif){
  
}



