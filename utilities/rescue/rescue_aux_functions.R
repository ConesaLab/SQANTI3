# function to search reference transcripts to see if they meet rescue criteria
rescue_lost_reference <- function(ref_id, classif){
  
  # filter classification
  classif_ref <- classif %>% 
    dplyr::filter(associated_transcript == ref_id)
  
  # check for FSM
  ref_check <- any(classif_ref$structural_category == "full-splice_match")
  
  # if there is an FSM associated to the lost reference, return lost reference
  if(ref_check == TRUE){
    ref_df <- tibble::enframe(ref_id, name = NULL, value = "isoform")
    return(ref_df)
    
  # if there is no FSM associated, return ISM artifacts
  }else{
    ism_df <- classif_ref %>%
      dplyr::select(isoform)
    
    return(ism_df)
  }
}


# function to create best_match_id column in rescue table
find_best_match_id <- function(candidate_id, rescue_table){
  
  # select only rows from rescue candidate
  rescue_sub <- rescue_table %>% 
    dplyr::filter(rescue_candidate == candidate_id)
    
  # select primary alignment
  sam <- rescue_sub %>% 
    dplyr::filter(sam_flag == 0)
      
    if(nrow(sam) == 1){
      id <- sam$mapping_hit
  
    }else{
      id <- paste(rescue_sub$mapping_hit, 
                  collapse = ",")
    }
  
  match_result <- tibble::tibble(rescue_candidate = candidate_id, 
                                 best_match_id = id)
  
  return(match_result)
}
