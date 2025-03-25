apply_rules <- function(isoform_info, force_multiexon){
  if (force_multiexon & as.numeric(isoform_info["exons"])==1){
  final_is_isoform=FALSE
  }else{
  sc = as.character(isoform_info["structural_category"])
  final_is_isoform=TRUE
  # detect if there is any specific rule for the SC of the isoform
  if (sc %in% names(json_df)) {
    final_is_isoform=FALSE
    # iterate all the independent rules for a certain SC
    for (p in which(names(rules_list)==sc)){
      is_isoform=TRUE
      rules <- rules_list[[p]]
      # iterate through the rules defined
      for (i in c(1:length(rules$rule))){
        parameter <- rules[i, "column"]
        if ( ! is.na(isoform_info[parameter])){ # if NA in the field, rule doesn't apply
          if (rules[i, "type"] == "Min_Threshold"){
            if (as.numeric(isoform_info[parameter]) < as.numeric(rules[i, "rule"])){
              is_isoform=FALSE
              break
            }
          }else if (rules[i, "type"] == "Max_Threshold"){
            if (as.numeric(isoform_info[parameter]) > as.numeric(rules[i, "rule"])){
              is_isoform=FALSE
              break
            }
          }else if (rules[i, "type"] == "Category"){
            cat_rules <- rules[rules$column == parameter, ]
            if ( ! tolower(isoform_info[parameter]) %in% cat_rules[,"rule"]){
              is_isoform=FALSE
              break
            }
          }
        }else{
          is_isoform=FALSE
          break
        }
      }
      final_is_isoform=final_is_isoform | is_isoform
    }
 # the isoform has a SC different, if will be evaluated with the rules of "rest" (if any was defined)  
  }else if ("rest" %in% names(json_df)){
    final_is_isoform=FALSE
    for (p in which(names(rules_list)=="rest")){
      is_isoform=TRUE
      rules <- rules_list[[p]]
      for (i in c(1:length(rules$rule))){
        parameter <- rules[i, "column"]
        if (! is.na(isoform_info[parameter])){
          if (rules[i, "type"] == "Min_Threshold"){
            if (as.numeric(isoform_info[parameter]) < as.numeric(rules[i, "rule"])){
              is_isoform=FALSE
              break
            }
          }else if (rules[i, "type"] == "Max_Threshold"){
            if (as.numeric(isoform_info[parameter]) > as.numeric(rules[i, "rule"])){
              is_isoform=FALSE
              break
            }
          }else if (rules[i, "type"] == "Category"){
            cat_rules <- rules[rules$column == parameter, ]
            if ( ! tolower(isoform_info[parameter]) %in% cat_rules[,"rule"]){
              is_isoform=FALSE
              break
            }
          }
        }else{
         is_isoform=FALSE
         break
        }
      }
      final_is_isoform=final_is_isoform | is_isoform
    }
  }
  }
  if (final_is_isoform){
   return("Isoform")
  }else{
   return("Artifact")
  }
}

### This function will apply the same filtering process but it will return the reasons for filtering a transcript. Use only on artifacts
get_reasons <- function(isoform_info, force_multiexon){
  if (force_multiexon & as.numeric(isoform_info["exons"])==1){
  reasons <- c("Mono-exon")
  }else{
  sc = as.character(isoform_info["structural_category"])
  reasons <- c()
  # detect if there is any specific rule for the SC of the isoform
  if (sc %in% names(json_df)) {
    # iterate all the independent rules for a certain SC
    for (p in which(names(rules_list)==sc)){
      rules <- rules_list[[p]]
      # iterate through the rules defined
      for (i in c(1:length(rules$rule))){
        parameter <- rules[i, "column"]
        if ( ! is.na(isoform_info[parameter])){ # if NA in the field, rule doesn't apply
          if (rules[i, "type"] == "Min_Threshold"){
            if (as.numeric(isoform_info[parameter]) < as.numeric(rules[i, "rule"])){
              reason_line <- paste("Low", parameter)
              reasons <- c(reasons, reason_line)
            }
          }else if (rules[i, "type"] == "Max_Threshold"){
            if (as.numeric(isoform_info[parameter]) > as.numeric(rules[i, "rule"])){
              reason_line <- paste("High", parameter)
              reasons <- c(reasons, reason_line)
              }
          }else if (rules[i, "type"] == "Category"){
            cat_rules <- rules[rules$column == parameter, ]
            if ( ! tolower(isoform_info[parameter]) %in% cat_rules[,"rule"]){
              reason_line <- paste("Out", parameter)
              reasons <- c(reasons, reason_line)
              }
          }
        } else {
          reason_line <- paste("NA value in", parameter)
          reasons <- c(reasons, reason_line)
        }
      }
    }
    # the isoform has a SC different, if will be evaluated with the rules of "rest" (if any was defined)  
  }else if ("rest" %in% names(json_df)){
    for (p in which(names(rules_list)=="rest")){
      rules <- rules_list[[p]]
      for (i in c(1:length(rules$rule))){
        parameter <- rules[i, "column"]
        if (! is.na(isoform_info[parameter])){
          if (rules[i, "type"] == "Min_Threshold"){
            if (as.numeric(isoform_info[parameter]) < as.numeric(rules[i, "rule"])){
              reason_line <- paste("Low", parameter)
              reasons <- c(reasons, reason_line)
            }
          }else if (rules[i, "type"] == "Max_Threshold"){
            if (as.numeric(isoform_info[parameter]) > as.numeric(rules[i, "rule"])){
              reason_line <- paste("High", parameter)
              reasons <- c(reasons, reason_line)
            }
          }else if (rules[i, "type"] == "Category"){
            cat_rules <- rules[rules$column == parameter, ]
            if ( ! tolower(isoform_info[parameter]) %in% cat_rules[,"rule"]){
              reason_line <- paste("Out", parameter)
              reasons <- c(reasons, reason_line)
            }
          }
        } else {
          reason_line <- paste("NA value in", parameter)
          reasons <- c(reasons, reason_line)
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






