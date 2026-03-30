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

# function to extract transcript and gene IDs from GTF
# Replacement for BUSpaRse::tr2g_gtf that is not available on ARM64 Macs
get_tr2g_gtf <- function(gtf_file) {
  # Read GTF file
  # GTF has 9 columns, no header
  # Use readr::read_tsv for speed and convenience
  gtf <- readr::read_tsv(gtf_file, col_names = FALSE, comment = "#", show_col_types = FALSE)
  
  # Filter for transcripts
  # Column 3 is feature type
  gtf_transcripts <- dplyr::filter(gtf, X3 == "transcript")
  
  if (nrow(gtf_transcripts) == 0) {
    gtf_transcripts <- dplyr::filter(gtf, X3 == "exon")
  }
  
  # Extract attributes from column 9 (X9)
  attributes <- gtf_transcripts$X9
  
  # Extract gene_id and transcript_id using regex
  # Assuming standard GTF format: gene_id "ID"; transcript_id "ID";
  
  gene_ids <- stringr::str_extract(attributes, 'gene_id "[^"]+"')
  gene_ids <- stringr::str_remove(gene_ids, 'gene_id "')
  gene_ids <- stringr::str_remove(gene_ids, '"')
  
  transcript_ids <- stringr::str_extract(attributes, 'transcript_id "[^"]+"')
  transcript_ids <- stringr::str_remove(transcript_ids, 'transcript_id "')
  transcript_ids <- stringr::str_remove(transcript_ids, '"')
  
  res <- data.frame(transcript = transcript_ids, gene = gene_ids, stringsAsFactors = FALSE)
  res <- dplyr::distinct(res)
  
  return(res)
}
