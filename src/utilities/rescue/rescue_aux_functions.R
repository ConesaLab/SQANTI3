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