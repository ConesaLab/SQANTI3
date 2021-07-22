##### Summary of ML variables in artifacts and isoforms #####

compare_MLvariables <- function(classification, variable, importance){
  
  require(ggplot2)
  require(magrittr)
  
  # Select variable for evaluation plot
  var_df <- classification %>% 
    dplyr::select(structural_category, filter_result,
                  dplyr::all_of(variable))
  
  # Get variable column info (class, name)
  var_type <- purrr::map_chr(var_df, class)
  var_name <- names(var_type[3])
  
  # Rename variable column to handle during plotting
  var_df <- var_df %>% dplyr::rename(variable = variable)
  
  # Round importance
  importance <- round(importance, 2)
  
  # Generate plot by type
  if(var_type[3] == "numeric"){
    
    p <- ggplot(var_df) +
      ggtitle(paste(var_name, "-", "ML importance:", importance)) +
      geom_boxplot(aes(x = structural_category, y = log(abs(variable)+1),
                       fill = filter_result), outlier.size = 0.2) +
      labs(x = "Structural category", y = paste0("log( |", var_name, "| +1)")) +
      theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
      scale_fill_manual("Filter result", 
                        values = c("deepskyblue3", "darkgoldenrod2"))
    
    return(p)
    
  } else{
    
    var_df <- var_df %>% 
      dplyr::filter(!is.na(variable)) %>% 
      dplyr::mutate(variable = factor(variable))
    
    p <- ggplot(var_df) +
      ggtitle(paste(var_name, "-", "ML importance:", importance)) +
      geom_bar(aes(x = filter_result, fill = variable), stat = "count", 
               width = 0.8, color = "black", position = "dodge") +
      labs(x = "Filter result", y = "Transcript no.") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
      RColorConesa::scale_fill_conesa(paste0(var_name), reverse = FALSE,
                                      palette = "complete", discrete = TRUE) +
      facet_grid(~structural_category, scales = "free")
    
    return(p)
    
  }
}
