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





