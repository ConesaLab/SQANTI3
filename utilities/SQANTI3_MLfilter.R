#!/bin/bash Rscript
#
##########################################################
#########     SQANTI3 MACHINE LEARNING FILTER    #########
##########################################################
#
# Authors: Cecile Pereira, Lorena de la Fuente, Francisco Pardo, Leandro Balzano-Nogueira, Ana Conesa
# Contact: f.pardo.palacios@gmail.com
# Affiliation: Institute for Integrative Systems Biology, CSIC, Valencia, Spain
#
# Last updated: June/30/2021
#
# Changes with respect to previous version
# RM instead of FSM are used as positive set
# New SQANTI features are used in the ML, related to CAGE, polyA and Ratio_TSS
# New options added
#   1. --output_directory



### Load required packages
library(caret)           # confusion matrix: caret_6.0-76 Version
library(dplyr) 
library(ggplot2)
library(ggpubr)
library(ipred)           # Improved Predictors package
library(MLmetrics)
library(nnet)            # neural networks
library(optparse)
library(partykit)        # "Elegant" tree design package
library(pROC)
library(randomForest)    # randomForest
library(ROCR)            # ROC curve
library(ROSE)
library(rpart)           # Rpart package
library(rpart.plot)      # Tree design 



### Definition of inputs

option_list = list(
  make_option(c("-c","--sqanti_classif"),type="character", default = NULL, 
              help="SQANTI classification output file"),
  make_option(c("-o","--output_directory"),type="character", default="SQANTI3_filter_out", 
              help="Output directory name"),
  make_option(c('-w','--wilcox'),type="integer",default = 0, 
              help="Feature selection option: 0 = max rfe value, 1 = minimum 
              number of feature with a mean non different to the mean of the 
              classifier obtain on all the features."),
  make_option(c("-t","--percent_training"),type="double",default = 0.8,
              help="the percentage of data that goes to training (parameter p of 
              the function createDataPartition)"),
  make_option(c("-p","--TP"), type="character",default = NULL,
              help="file containing the list of the TP transcripts, 
              one ID by line, no header"),
  make_option(c("-n","--TN"), type="character",default = NULL,
              help="file containing the list of the TN transcripts, 
              one ID by line, no header"),
  make_option(c("-j","--threshold"), default = "0.7",
              help="machine learning probability threshold to classify 
              positive isoforms"),
  make_option(c("-i","--intrapriming"), default = "80", 
              help="polyA percentage thereshold to flag an isoform as 
              intra-priming"),
  make_option(c("-f","--force_fsm_in"), default = FALSE, 
              help="FMS transcripts are not filtered")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser) #list of the args
opt$threshold = as.numeric(opt$threshold)
opt$intrapriming = as.numeric(opt$intrapriming)
dir.create(opt$output_directory)

# Reading data
d <- read.table(file = opt$sqanti_classif, sep ="\t", header= TRUE,as.is =TRUE)

# Isoform name column goes to rownames and is deleted from the table
rownames(d) = d$isoform
d = d[,-which(colnames(d) == "isoform")]

# Check data for mono-exons. Only multi-exon transcripts are subject to ML filter
monoexons <- nrow(d[which(d$exons == 1),])  # number of monoexon transcripts
multiexons <- nrow(d) - monoexons

if(multiexons == 0){
  message("Machine Learning filtering won't be applied because all the isoforms are mono-exon")
  flag = FALSE
}else{
  message("Number of multi-exon transcript isoforms:")
  message(multiexons)
}

# Create the set of TP and TN
if(is.null(opt$TP) | is.null(opt$TN)){
  message("Training set not provided, it will be created from input data.")
  FSM_set <-  rownames(d[d$structural_category == "full-splice_match" & d$exons > 1,])
  RM_set <-   rownames(d[d$subcategory == "reference_match",])
  NNIC.NC_set <- rownames(d[(d$structural_category == "novel_not_in_catalog" & 
                               d$all_canonical == "non_canonical"),])
  flag = TRUE
  if (nrow (NNIC.NC_set) > 40 ) {
    Negative_set <- NNIC.NC_set
    if (nrow(FSM_set) > 40 ) {
      Positive_set <- FSM_set
      if (nrow (RM_set) > 40) {
        Positive_set <- RM_set
      } else { message ("Not enough Reference Match transcritps, 
                      Full-Splice-Match transcripts will be used as Positive set.")}
    } else { message ("Not enough Full-Splice-Match transcripts, ML filter is skipped.")
            flag = FALSE
    }
  } else { 
    message  ("Not enough Novel_Not_in_Catalog + Non_Cannonical transcripts, 
            ML filter is skipped.")
  flag = FALSE 
  }
}


#############################
###### Machine Learning #####
#############################
 
if (flag) {
  message("Initializing Machine Learning computations")
  
### Preparation of classification table for ML 

  # Create new classification table for ML data preparation
  d1 <- d
  
  # Sum all the columns that contain FL reads and create only one column
  if (any(grep('^FL\\.', names(d1)))){d1$FL= rowSums(d1[,grep('^FL\\.', names(d1))])}
  
  # Case special of polyA motif. When present, simply replace by "Motif_found"
  d1[which(!is.na(d1$polyA_motif)), "polyA_motif"] <- "Motif_found"
  
  #  Change NA in columns by an appropriate replacement
  NA_columns <- c("polyA_motif", "within_cage_peak", 'n_indels', "n_indels_junc", 
                  "predicted_NMD", "min_sample_cov", "min_cov", "ratio_exp", "bite", 
                  "diff_to_gene_TSS", "diff_to_gene_TTS" , "dist_to_polya_site", 
                  "dist_to_cage_peak", 'within_polya_site', "polyA_dist")
  replacement.na <- c("No_polyA_motif", 0,0, 0, "non_coding",0, 0,0, FALSE, 
                      -11000, -11000, -11000, -11000, FALSE, -11000)
  for (i in 1: length(NA_columns)) {
    sel.column <- which(colnames (d1) == NA_columns [i])
    d1[which(is.na(d1[,sel.column])), sel.column] <- replacement.na[i]
  }
  
 
  # Case special sdCov, If all NAs, replace by 0. If some NA replace with 
  # median value of sdCov in the dataset
  if (all(is.na(d1$sd_cov))){
    d1$sd_cov = 2
  } else{
    medt2 = median(as.numeric(d1[!is.na(d1$sd_cov), "sd_cov"]))
    d1[is.na(d1$sd_cov),"sd_cov"] <- medt2
  }
  
  # Convert in factors columns with categorical variables
  categorical <- c("FSM_class", "coding", "bite", "within_cage_peak", 
                   "polyA_motif" , "within_polya_site", "predicted_NMD")
  for (x in categorical){
    d1[,x] <- as.factor(d1[,x])
  }
  
  # Convert in integers columns with numerical variables
  integers <- c("diff_to_gene_TSS", "diff_to_gene_TTS", "min_sample_cov", 
                "min_cov", "ratio_exp" , "polyA_dist", "dist_to_cage_peak")
  for (x in integers){
    d1[,x] <- as.integer(d1[,x])
  }
  
  r = as.vector(which(apply(d1,2,function(x) (anyNA(x)))))
  if (length(r)>0){  d1 = d1[,-r] }
  
  # Removing columns with zero variance
  message("-------------------------------------------------")
  message("Preprocessing: removing variables with near zero variance")
  nzv = nearZeroVar(d1)
  if(length(colnames(d1)[nzv]) != 0){
    message("removed columns: ")
    message(colnames(d1)[nzv])
    d1 = d1[,-nzv]
  } else {
    message("Preprocessing: no variables with near-zero variance found")}
  
  # Calculating redundant variables
  message("Preprocessing: removing the features with a correlation > 0.9")
  r = as.vector(which(apply(d1,2,function(x) (anyNA(x)))))
  if (length(r) > 0){d1 = d1[,-r] }
  d2 <- d1[,sapply(d1,class)%in%c("numeric", "integer")] # selecting only numeric variables
  descrCorr = cor(d2)
  highCorr <- findCorrelation(descrCorr, cutoff = 0.9, verbose = TRUE, names=TRUE)
  if(length(highCorr)>0){
    message("List of removed features: ")
    message ("highCorr")
    d1 = d1[,-which(colnames(d1)%in%highCorr)]
  } else {
    message("No feature removed")
  }
  
  # END  preparing data for ML 
  ############################
  
  
  ##########################
  ####  ML algorithm   #####
  ##########################
  
  #### Creating postive and negative set
  message("-------------------------------------------------")
  message("Creating positive and negative sets for training and testing")
  
  
  # Removing monoexons
  dmult <- d1[which(d1$exons != 1),] #considering only multiexon 
  
  
  #### Making the training set
  
   if(is.null(opt$TP)){
    Positive_df <- dmult[intersect(rownames(dmult), rownames(Positive_set)),colnames(dmult)]
    Negative_df <- dmult[intersect(rownames(dmult), rownames(Negative_set)),colnames(dmult)]
    trainingset = rbind(Positive_df,Negative_df)
    Class = factor(c(rep("POS",nrow(Positive_set)),rep("NEG",nrow(Negative_set))))
  } else {
    TP = read.table(opt$TP,as.is = TRUE)
    TN = read.table(opt$TN,as.is = TRUE)
    Class = factor(c(rep("POS",length(TP$V1)),rep("NEG",length(TN$V1))))
    trainingset = dmult[intersect(rownames(dmult), c(TP$V1,TN$V1)),]
  }
  
  # Remove columns that are not informative for ML
  colRem = c('chrom','strand','associated_gene', 'associated_transcript', 
             'ref_length','ref_exons', 'ORF_length', 'CDS_length', 'CDS_start', 
             'CDS_end', 'CDS_genomic_start', 'CDS_genomic_end', 'all_canonical',
             'seq_A_downstream_TTS', 'ORF_seq', 'subcategory', 'structural_category')
  trainingset = trainingset[,-which(colnames(trainingset)%in%colRem)]
  
  
  ####  Partition
  message("Partition train set / test set")
  message("Proportion of the label data used for the training (between 0 and 1):")
  message(opt$percent_training)
  
  set.seed(123)
  inTraining = createDataPartition(Class,p = opt$percent_training,list = FALSE, times = 1)[,1]
  
  training = trainingset[inTraining,]
  testing = trainingset[-inTraining,]
  
  message("Description of the training set:")
  message("Table number of positive and negative examples in the training set")
  message(table(d[d[,1]%in%rownames(training),]$structural_category))
  message("Table number of positive and negative examples in the testing set")
  message(table(d[d[,1]%in%rownames(testing),]$structural_category))
  
  #10 times 10 cross validation
  ctrl = trainControl(method = "repeatedcv", repeats = 10,
                    classProbs = TRUE,
                    sampling = 'down',returnData = TRUE,
                    savePredictions = TRUE, returnResamp ='all')
  message("-------------------------------------------------")
  message("Training Random Forest Classifier. This can take up to several hours")
  message("Random Forest parameters:")
  message("-Down sampling in training set")
  message("-10 cross validation")
  
  set.seed(1)
  
  #default: 500 trees
  randomforest <- caret::train(x = training, y = as.factor(Class[inTraining]),
                               method ='rf',
                               tuneLength = 15,
                               metric = "Accuracy",
                               trControl = ctrl)
  
  save(randomforest, file = "randomforest.RData")
  load("randomforest.RData")
  

  #### Apply classifier on test set

  message("-------------------------------------------------")
  message("Classifier performance in test set")
  test_pred_prob = predict(randomforest,testing,type = 'prob')
  pred = factor(ifelse(test_pred_prob$POS >= opt$threshold,"POS","NEG"))
  a = data.frame(POS=test_pred_prob$POS,
                 NEG=test_pred_prob$NEG,
                 pred, 
                 obs = Class[-inTraining])
  rownames (a) <- rownames(test_pred_prob)
  
  message("AUC, Sens, Spec on the test set")
  message(twoClassSummary(a,lev = levels(a$obs)))
  write("AUC, Sens, Spec on the test set",
        file = paste(opt$output_directory,'/statistics_testSet.txt',sep = ''))
  write(twoClassSummary(a,lev=levels(a$obs)),
        file = paste(opt$output_directory,'/statistics_testSet.txt',sep=''), append = TRUE)
  write.table(a,paste(opt$output_directory,'/Pred_test_and_class.txt',sep =''), 
              quote = F, sep = "\t", row.names = TRUE)
  
  cm = confusionMatrix(data = pred,
                       reference = Class[-inTraining],
                       positive = "POS")
  info = "rows:predictions \ncolumns:reference"
  write.table(paste(info,"\n"),
              file = paste(opt$output_directory,"/confusionMatrix_testSet.txt",sep=''), 
              quote = F, col.names = F, row.names = F)
  write.table(cm$table,
              file = paste(opt$output_directory,"/confusionMatrix_testSet.txt",sep=''), 
              append = TRUE, quote = F, col.names = T, row.names = T)
  write.table("\n",file = paste(opt$output_directory,"/confusionMatrix_testSet.txt",sep=''),  
              append = TRUE,quote = F, col.names = F, row.names = F)
  write.table(cm$overall,file = paste(opt$output_directory,"/confusionMatrix_testSet.txt",sep=''),
              append = TRUE, quote = F, col.names = F)
  write.table(cm$byClass,file = paste(opt$output_directory,"/confusionMatrix_testSet.txt",sep=''),
              append = TRUE, quote = F, col.names = F)
  
  # Varaible importance for the prediction
  imp = varImp(randomforest,scale = FALSE)
  imp <- imp$importance
  imp <- data.frame(rownames(imp),imp) 
  imp <- imp[order(-imp$Overall, decreasing = FALSE),]
  imp <- imp[,-1, drop = FALSE] # added Leo
  write.table(imp,file = paste(opt$output_directory,
                               "/VariableImportanceInClassifier.txt",sep =''), 
              sep = "\t", quote = F, col.names = F)
  par(mar = c(5.1, 10.1, 4.1, 2.1)) # all sides have 3 lines of space
  pdf("variable.importance.pdf")
  barplot (as.matrix(t(imp)) , horiz = TRUE, las = 2, col = "lightblue", 
           main = "Variable Importance In Classifier")
  dev.off()
  
  ###### ROC curve
  # 1) in function of the probability on the test set 
  # (not the same proportion of positives and negatives)
  fileroc = paste(opt$output_directory,"/ROC_curve_testset.pdf",sep='')
  pdf(fileroc)
  r = roc(as.numeric(Class[-inTraining]),test_pred_prob$POS,percent = TRUE)
  auc(r)
  plot.roc(r, main = "ROC with unbalanced classes")
  text(20,10,paste('AUC =',signif(auc(r),4)))
  text(20,5,paste('CI 95% = [',signif(ci(r)[1],4),',',signif(ci(r)[2]),']'))
  
  # # 2) same proportion positives and negatives on the test set:
  #testing
  #list of the testing positives
  alltestpos = which(Class[-inTraining] == 'POS')
  alltestneg = which(Class[-inTraining] == 'NEG')
  nbpos = length(alltestpos)
  nbneg = length(alltestneg)
  if(nbpos < nbneg){
    sampleneg = sample(alltestpos,nbpos,replace = FALSE)
    newtest = testing[c(sampleneg,alltestpos),]
    Classnewtest = factor(c(rep('NEG',nbpos),rep('POS',nbpos)))
  } else {
    samplepos = sample(alltestpos,nbneg,replace=FALSE)
    newtest = testing[c(samplepos,alltestneg),]
    Classnewtest = factor(c(rep('POS',nbneg),rep('NEG',nbneg)))
  }
  set.seed(1)
  test_pred_prob2 = predict(randomforest, newtest, type = 'prob')
  
  r = roc(as.numeric(Classnewtest),test_pred_prob2$POS,percent = TRUE)
  auc(r)
  #Area under the curve: 99.29%
  plot.roc(r, main = "ROC with balanced classes")
  lines(r,col = 'red')
  ci(r)
  text(20,20,paste('AUC =',signif(auc(r),4)),col = 'red')
  text(20,15,paste('CI 95% = [',signif(ci(r)[1],4),',',signif(ci(r)[2]),']'),col = 'red')
  dev.off()
  
  ##################################
  #### Applying classifier to data #
  ##################################
  message("-------------------------------------------------")
  message("Classifier in our dataset")
  
  isoform.predict = predict(randomforest,dmult[,colnames(training)],type = 'prob')
  colnames(isoform.predict) = gsub("NEG","NEG_MLprob", colnames(isoform.predict))
  colnames(isoform.predict) = gsub("POS","POS_MLprob", colnames(isoform.predict))
  message( "Random forest prediction done")
  
  ## Adding predictions to classification table
  classified.isoforms = cbind(dmult[rownames(isoform.predict),],isoform.predict)
  
  negatives = classified.isoforms[classified.isoforms$POS_MLprob < opt$threshold,]
  if (opt$force_fsm_in) {
    negatives <- negatives[negatives$structural_category != "full-splice_match",]
  }
  
  classified.isoforms$"ML_classifier" <- "Positive"
  classified.isoforms$"ML_classifier"[rownames(classified.isoforms) %in% 
                                        rownames(negatives)] <- "Negative"
  message( "Results")
  message(table(classified.isoforms$"ML_classifier"))
  
  ####################################
  ## Add monoexons back to dataset  ##
  ####################################
  dme <- d1[setdiff(rownames(d1), rownames(isoform.predict)),]
  dme$"NEG_MLprob"  <- dme$"POS_MLprob" <- dme$"ML_classifier"  <- NA
  
  d1 <- rbind (classified.isoforms, dme[,colnames(classified.isoforms)])
}


####################################
###### INTRA-PRIMING FILTERING #####
####################################

message("-------------------------------------------------")
message("Intrapriming filtering in our dataset")

d1[,"intra_priming"] = d1$perc_A_downstream_TTS > as.numeric(opt$intrapriming) & 
  !(d1$structural_category %in% c("full-splice_match") )

message("Intra-priming filtered transcripts:")
message(table(d1$'intra-priming'))

d1  <- d1[rownames(d),] # reordering of d and d1 to have transcripts in the same order

##### Final result tables
##########################


# d: Initial table with new column that indicates if transcript is isoform or artifact
d$SQANTI_filter <- "Isoform"
d[which(d1$ML_classifier == "Negative" | 
          d1$intra_priming == TRUE),"SQANTI_filter"] <- "Artifact"
write.table(d,file = paste(opt$output_directory,
                           "/SQANTI_classification_ML_prediction.txt",
                           sep =''),
            quote = FALSE, col.names = TRUE, sep ='\t',row.names = TRUE)

# d1: Table with variables modified by the ML filter and with the result of ML 
# filter and intra-priming evaluation
d1 = data.frame(d1, SQANTI_filter = d$SQANTI_filter)
write.table(d1,file = paste(opt$output_directory,
                            "/Output_ML_filter_agorithm.txt",
                            sep = ''),
            quote = FALSE, col.names = TRUE, sep ='\t',row.names = TRUE)


# new_classification: classification table with only Isoform transcripts
new_classification <- d[d$"SQANTI_filter" == "Isoform",]
write.table(new_classification, 
            file = paste(opt$output_directory,
                         "/SQANTI_classification_ML_filtered.txt",
                         sep =''),
            quote = FALSE, col.names = TRUE, sep ='\t',row.names = TRUE)

# discarded_transcripts and why
discarded <- d1[which(d1$ML_classifier == "NEG" | d1$`intra-priming`== TRUE),]
write.table(discarded, file = paste(opt$output_directory,
                                    "/Discarded_transcripts.txt",
                                    sep=''),
            quote = FALSE, col.names = TRUE, sep ='\t',row.names = TRUE)
    
message("-------------------------------------------------")
message("*** SQANTI3 ML filter finished!")

################################################



##### Summary of results #####
# plot.compare is a function to evaluate the values of the features used in the 
# ML filter in Artifacts and Isoforms per structural category

i = 2
classification = d1; myfeature = features.eval[i]; imp = importance[i]


plot.compare <- function (classification, myfeature, imp) {
  library(reshape2)
  mycategories = c("full-splice_match", "incomplete-splice_match", "novel_in_catalog",
                   "novel_not_in_catalog", "intergenic", "fusion", "genic", 
                   "antisense", "genic_intron")
  mylabels = c("FSM", "ISM", "NIC", "NNC", "Inter", "Fus", "Gen", "Anti", "Intron")
  col <- which(colnames(classification) == myfeature)
  
  
  if (class(classification[,col]) != "factor") {
    if (myfeature != "count")  {
      dat <- data.frame(category = classification$structural_category, 
                        filter = classification$SQANTI_filter, 
                        feature = classification[,col])
      dat2 <- with(dat, tapply(feature,list(category,filter), median))
    } else {
      dat <- data.frame(category = classification$structural_category, 
                        filter = classification$SQANTI_filter)
      dat2 <- table(dat)
      dat2 <- dat2[mycategories,]
    }
    dat3 <- melt(dat2) ; names(dat3) <- c("category", "filter", "feature")
    dat3$category <- factor(dat3$category, levels = mycategories)
    dat3$filter <- factor(dat3$filter, levels = c("Isoform", "Artifact"))
    
    p <- ggplot(dat3, aes(x = category, y =  feature)) +
      geom_bar(
        aes(color = filter, fill = filter),
        stat = "identity", position = position_dodge(0.8),
      ) +
      scale_color_manual(values = c("#0073C2FF", "#EFC000FF")) +
      scale_fill_manual(values = c("#0073C2FF", "#EFC000FF")) +
      ggtitle(paste(myfeature, "importance", imp, sep = " ")) + 
      scale_x_discrete(breaks = mycategories,
                       labels = mylabels)
    
  } else {   # code for categorical features
    col <- which(colnames(classification) == myfeature)
    dat <- data.frame(category = classification$structural_category, 
                      filter = classification$SQANTI_filter, feature = classification[,col])
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
  
  print(p)
}

features.eval <- c("count", rownames(imp), "intra_priming")
importance <- c("NA", round(imp[,1],2), "NA")
pdf("Evaluation.ML.pdf")
for ( i in 1:length(features.eval)) {
  print(features.eval[i])
  plot.compare (classification = d1, myfeature = features.eval[i], imp = importance[i])
}
dev.off()


plot(density(log(classification$ratio_TSS[classification$structural_category == "incomplete-splice_match"])), col = "red")
lines(density(log(classification$ratio_TSS[classification$structural_category == "full-splice_match"])))
lines(density(log(classification$ratio_TSS[(classification$structural_category == "incomplete-splice_match" & classification$SQANTI_filter == "Isoform")])), col = "blue")





