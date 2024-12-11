#!/bin/bash Rscript
#
##########################################################
#########     SQANTI3 MACHINE LEARNING FILTER    #########
##########################################################
#
# Authors: Cecile Pereira, Lorena de la Fuente, Francisco Pardo, 
# Leandro Balzano-Nogueira, Ana Conesa, Ángeles Arzalluz-Luque
# Contact: f.pardo.palacios@gmail.com
# Affiliation: Institute for Integrative Systems Biology, CSIC, Valencia, Spain
#
# Last updated: July/13/2021
#
# Changes with respect to previous version
# - RM instead of FSM are used as positive set.
# - New SQANTI features are used in the ML, i.e. related to CAGE,  and Ratio_TSS.
# - Outputs: 
#   1. Classification table + ML results.
#   2. Single-column text file with IDs of transcripts labeled as isoforms by
#     ML filter and intra-priming prediction.
# - New options added:
#   1. --output and --dir to control file prefix and output directory, respectively.
#   2. --intermediate_files to output ML-modified classification table.
#   3. NNC as true negative set when no NNC non-canonical transcripts are available.



### Make script argument list

option_list = list(
  optparse::make_option(c("-c","--sqanti_classif"), type="character", default = NULL, 
              help="SQANTI classification output file."),
  optparse::make_option(c("-o","--output"), type="character", default = "SQANTI3", 
              help="Output classification file prefix."),
  optparse::make_option(c("-d","--dir"), type="character", 
              help="Output directory."),
  optparse::make_option(c("-t","--percent_training"),type="double",default = 0.8,
              help="Default: 0.8. Proportion of the data that goes to training 
              (parameter p of the function createDataPartition)"),
  optparse::make_option(c("-p","--TP"), type="character",default = NULL,
              help="Path to file containing the list of the TP transcripts, 
              one ID by line, no header (optional). If not supplied, it will be 
              generated from input data."),
  optparse::make_option(c("-n","--TN"), type="character",default = NULL,
              help="Path to file containing the list of the TN transcripts, 
              one ID by line, no header (optional). If not supplied, it will be 
              generated from input data."),
  optparse::make_option(c("-j","--threshold"), type="double", default=0.7,
              help="Default: 0.7. Machine learning probability threshold to classify, 
              transcripts as positive isoforms."),
  optparse::make_option(c("-i","--intrapriming"), type="integer", default=60, 
              help="Default: 60.  thereshold (i.e. number of A's downstream 
              the TTS) to flag an isoform as an intra-priming artifact."),
  optparse::make_option(c("-f","--force_fsm_in"), type="logical", default = FALSE, 
              help="Default: FALSE. When TRUE, forces retaining FMS transcripts 
              regardless of ML filter, FSM are threfore not filtered."),
  optparse::make_option(c("-e", "--force_multi_exon"), type="logical", default = FALSE,
              help="Default: FALSE. When TRUE, forces retaining only multi-exon 
              transcripts, all mono-exon isoforms will be automatically removed."),
  optparse::make_option(c("-m", "--intermediate_files"), type="logical", default=FALSE,
              help="Default: FALSE. When TRUE, outputs ML filter intermediate 
              files."),
  optparse::make_option(c("-r", "--remove_columns"), type="character", default=NULL,
              help="Path to single-column file (no header) containing the names of 
              the columns in SQ3's classification.txt file that are to be excluded 
              during random forest training (optional)."),
  optparse::make_option(c("-z", "--max_class_size"), type="numeric", default=3000,
              help="Default: 3000. Maximum number of isoforms to include in True 
              Positive and True Negative sets. TP and TN sets will be downsized 
              to this value if they are larger.")
)

# Parse and handle provided arguments
opt_parser = optparse::OptionParser(option_list = option_list)
opt = optparse::parse_args(opt_parser) #list of the args
opt$threshold = as.numeric(opt$threshold)
opt$intrapriming = as.numeric(opt$intrapriming)



###############################
###### Generate ML inputs #####
###############################

message("-------------------------------------------------")
message("\n \t SQANTI3 Machine Learning filter")
message("\n--------------------------------------------------")

### Print specified parameters
message("\nCURRENT ML FILTER PARAMETERS:\n")
print(paste0(names(opt), ": ", opt))


### Read input data
message("\n\tINITIAL ML CHECKS:")
message("\nReading SQANTI3 *_classification.txt file...")

d <- read.table(file = opt$sqanti_classif, sep ="\t", header= TRUE,as.is =TRUE)

# Isoform name column goes to rownames and is deleted from the table
rownames(d) = d$isoform
d = d[,-which(colnames(d) == "isoform")]



# Check data for mono-exons. Only multi-exon transcripts are subject to ML filter
message("\nChecking data for mono and multi-exon transcripts...")

monoexons <- nrow(d[which(d$exons == 1),])  # number of monoexon transcripts
multiexons <- nrow(d) - monoexons

  if(multiexons == 0){
    message("\n\tWarning message: \n \t All isoforms in SQ3 classification file are mono-exon: 
            skipping ML filter.")
    
    run_ML = FALSE
    
  } else{
    message("\n \t ***Note: ML filter can only be applied to multi-exon transcripts. ")
    message(paste("\n \t", multiexons, "multi-exon transcript isoforms found in SQ3 classification file."))
  
  }



#### Check inputs for TP and TN sets

message("\nChecking input data for True Positive (TP) and True Negative (TN) sets...")

# Positive and negative sets not defined yet
Positive_set <- NULL
Negative_set <- NULL

# First check whether TP and TN sets are supplied (first if)
# If not available, create the set of TP and TN from input data
if(is.null(opt$TP) == FALSE & is.null(opt$TN) == FALSE){
  
  message("\n--TP and --TN arguments provided: using supplied set of isoforms as training set.")
  
  TP <- read.table(opt$TP)
  TN <- read.table(opt$TN)
  
  # keep only multi-exon transcripts in user-defined TP and TN sets
  tp_keep <- d[unlist(TP),]$exons > 1
  TP <- TP[tp_keep, , drop = FALSE]
  
  tn_keep <- d[unlist(TN),]$exons > 1
  TN <- TN[tn_keep, , drop = FALSE]
  
  
  if(nrow(TN) >= 250 & nrow(TP) >= 250){
    
    run_ML <- TRUE
    
    # define TP
    Positive_set <- unname(unlist(TP))
    message(paste0("\n\t - Total isoforms in user-defined TP set: ", 
                   length(Positive_set)))
    
    # define TN
    Negative_set <- unname(unlist(TN))
    message(paste0("\n\t - Total isoforms in user-defined TN set: ", 
                   length(Negative_set)))
    
  } else{
    
    run_ML <- FALSE
    
    message("\nWarning message:
            user-defined TP and TN sets must have >=250 isoforms! Skipping ML filter.")
  }
  
  
} else{
  message("\n\tWarning message: \n \t Training set not provided -will be created from input data.")
  
  FSM_set <- rownames(d[d$structural_category == "full-splice_match" & d$exons > 1,])
  RM_set <- rownames(d[d$subcategory == "reference_match" & d$exons > 1,])
  NNC_set <- rownames(d[d$structural_category == "novel_not_in_catalog",])
  NNC.NC_set <- rownames(d[(d$structural_category == "novel_not_in_catalog" & 
                               d$all_canonical == "non_canonical"),])
  run_ML = TRUE
  
  
  # 1. CHECK NEGATIVE SET REQUIREMENTS
  # Check whether number of NNC non-canonical is sufficient to run ML filter
  # If not, check NNC and see if it meets length requirement
  if (length (NNC.NC_set) < 250) {
    
    message("\nWarning message:
            \nNot enough (< 250) Novel Not in Catalog (NNC) + non-canonical transcripts.")

    if(length(NNC_set) >= 250){
      Negative_set <- NNC_set
      
      message("\nUsing Novel Not in Catalog (NNC) transcripts as True Negatives for training.")
      message(paste("\n \t - Total NNC isoforms: ", length(Negative_set)))
      
    } else{
      run_ML = FALSE 
      
      message("\nWarning message:
            \nNot enough (< 250) Novel Not in Catalog (NNC) transcripts, skipping ML filter.")
      message("\n\t***Note: try re-running ML filter with a user-defined TN set >=250 isoforms!")
      
    }
    
  } else { 
    Negative_set <- NNC.NC_set
    
    message("\nUsing Novel Not In Catalog non-canonical isoforms as True Negatives for training.")
    message(paste("\n \t - Total NNC non-canonical isoforms:", length(Negative_set)))
    
  }
  
  # 2. IF NEGATIVE_SET HAS BEEN DEFINED, RUN CHECKS TO CREATE POSITIVE_SET
  # Check whether RM set meets length requirement
  # If not, check FSM set for length and assign as TP
  if(is.null(Negative_set) == FALSE){
    
      if(length (RM_set) >= 250){
        Positive_set <- RM_set
        
        message("\nUsing FSM Reference Match isoforms as True Positives for training")
        message(paste("\n \t - Total reference match isoforms (FSM subcategory):", 
                      length(Positive_set)))
        
        # If RM are TP set, exclude diff_to_* columns
        colRem_RM <- c("diff_to_gene_TSS", "diff_to_gene_TTS",
                       "diff_to_TSS", "diff_to_TTS")
        
        message("\nExcluding diff_to_* columns before model training to prevent overfitting!")
        message(paste("\n \t The following columns will be added to the column removal list:", 
                      colRem_RM))
        
      }
      else if(length(FSM_set) >= 250){ 
        Positive_set <- FSM_set
        
        message("\nNot enough (< 250) Reference Match transcript isoforms among FSM, 
                      all FSM transcripts will be used as Positive set.")
        message(paste("\n \t - Total FSM isoforms:", length(Positive_set)))
        
      }
      else{ 
        message ("Warning message: 
                         \nNot enough (< 250) Full-Splice-Match transcripts, skipping ML filter.")
        message("\n\t***Note: try re-running ML filter with a user-defined TP set >=250 isoforms!")
        
        run_ML = FALSE
      } 
  }
}


#### Balance isoform number in TP and TN sets
if(length(Positive_set) != length(Negative_set)){
  
  message("\n\nBalancing number of isoforms in TP and TN sets...")
  
  # select minimum size among TP and TN sets
  sub_size <- min(length(Positive_set), length(Negative_set))
  
  message(paste0("\n\tMinimum set size: ", sub_size, " transcripts."))
  
      #### Downsize TP and TN sets if isoform no. is too large
      if(sub_size > opt$max_class_size){
        
        message(paste0("\nWarning message: min. set size is larger than --max_class_size!"))
        message(paste0("\n[!] TP and TN sets will be downsized to ", opt$max_class_size, 
                       " isoforms."))
        
        sub_size <- opt$max_class_size
      }
  
  message(paste0("\n\tSampled ", sub_size, " transcripts to define final TP and TN sets."))
  
  # sample sub_size number of isoforms as TP and TN
  Positive_set <- sample(Positive_set, sub_size)
  Negative_set <- sample(Negative_set, sub_size)
}


#### If not provided, output generated list of TP and TN
if(is.null(opt$TP) == TRUE & is.null(opt$TN) == TRUE){
  
  # create single-column tables
  Negative_df <- data.frame(transcripts = Negative_set)
  Positive_df <- data.frame(transcripts = Positive_set)
  
  # output tables
  write.table(Negative_df, 
              file = paste0(opt$dir, "/", opt$output, "_TN_list.txt"),
              quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)
  
  write.table(Positive_df, 
              file = paste0(opt$dir, "/", opt$output, "_TP_list.txt"),
              quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE)
  
  message("\nWrote generated TP and TN lists to files:")
  message(paste0("\n\t", opt$dir, "/", opt$output, "_TP_list.txt"))
  message(paste0("\n\t", opt$dir, "/", opt$output, "_TN_list.txt"))
}


#############################
###### Machine Learning #####
#############################
 
# Create new classification table for ML data preparation
d1 <- d


if (run_ML == TRUE) {
  message("\n-------------------------------------------------")
  message("\n\tML DATA PREPARATION:")
  
### Preparation of classification table for ML 
  
  # Sum all the columns that contain FL reads and create only one column
  message("\nAggregating FL counts across samples (if more than one sample is provided)...")
  
  if (any(grep('^FL\\.', names(d1)))){d1$FL= rowSums(d1[,grep('^FL\\.', names(d1))])}
  
  #  Change NA in columns by an appropriate replacement
  message("\nReplacing NAs with appropriate values for ML...")
  
  NA_columns <- c("within_CAGE_peak", 'n_indels', "n_indels_junc", 
                  "predicted_NMD", "min_sample_cov", "min_cov", "ratio_exp", "bite", 
                  "diff_to_gene_TSS", "diff_to_gene_TTS" , "dist_to_polyA_site", 
                  "dist_to_CAGE_peak", 'within_polyA_site', "polyA_dist",
                  "ratio_TSS")
  
  replacement.na <- c(0, 0, 0, "non_coding",0, 0, 0, FALSE, 
                      -11000, -11000, -11000, -11000, FALSE, -11000, 1)
  
  for (i in 1: length(NA_columns)) {
    sel.column <- which(colnames (d1) == NA_columns [i])
    d1[which(is.na(d1[,sel.column])), sel.column] <- replacement.na[i]
  }
  
 
  # Special case for sdCov: If all NAs, replace by 0. If there are NA values, 
  # replace with  median value of sdCov in the dataset
  if (all(is.na(d1$sd_cov))){
    d1$sd_cov = 2
  } else{
    medt2 = median(as.numeric(d1[!is.na(d1$sd_cov), "sd_cov"]))
    d1[is.na(d1$sd_cov),"sd_cov"] <- medt2
  }
  
  # Convert in factors columns with categorical variables
  message("\nHandling factor columns...")
  
  categorical <- c("FSM_class", "coding", "bite", "within_CAGE_peak", 
                   "polyA_motif_found" , "within_polyA_site", "predicted_NMD")
  for (x in categorical){
    d1[,x] <- as.factor(d1[,x])
  }
  
  # Convert in integers columns with numerical variables
  message("\nHandling integer columns...")
  
  integers <- c("diff_to_gene_TSS", "diff_to_gene_TTS", "min_sample_cov", 
                "min_cov", "ratio_exp" , "polyA_dist", "dist_to_CAGE_peak")
  for (x in integers){
    d1[,x] <- as.integer(d1[,x])
  }
  
  r = as.vector(which(apply(d1,2,function(x) (anyNA(x)))))
  if (length(r)>0){  d1 = d1[,-r] }
  
  # Removing columns with zero variance
  message("\nRemoving variables with near-zero variance...")
  nzv = caret::nearZeroVar(d1)
  if(length(colnames(d1)[nzv]) != 0){
    message("\tRemoved columns: ")
    print(paste(colnames(d1)[nzv]))
    d1 = d1[,-nzv]
  } else {
    message("\tNo variables with near-zero variance found.")}
  
  # Calculating redundant variables
  message(paste0("\nRemoving highly correlated features... (correlation threshold = 0.9).\n"))
  r = as.vector(which(apply(d1,2,function(x) (anyNA(x)))))
  if (length(r) > 0){d1 = d1[,-r] }
  d2 <- d1[,sapply(d1,class)%in%c("numeric", "integer")] # selecting only numeric variables
  descrCorr = cor(d2)
  highCorr <- caret::findCorrelation(descrCorr, cutoff = 0.9, verbose = TRUE, names=TRUE)
  message("\n\tList of removed features: ")
  if(length(highCorr)>0){
    message(paste("\t", highCorr))
    d1 = d1[,-which(colnames(d1) %in% highCorr)]
  } else {
    message("\tNo features removed.")
  }
  
  # END  of ML data preparation 
  ############################
  
  
  ##########################
  ####  ML algorithm   #####
  ##########################
  
  #### Creating positive and negative set
  message("\n-------------------------------------------------")
  message("\n\tRANDOM FOREST ALGORITHM RUN:")
  message("\nCreating positive and negative sets for classifier training and testing...")
  
  
  # Remove mono-exon transcripts from training set
  dmult <- d1[which(d1$exons != 1),]
  
  
  # Subset pre-processed data table (dmult) to make training dataset
  trainingset = rbind(dmult[Positive_set, ], 
                        dmult[Negative_set,])
  Class = factor(c(rep("POS", length(Positive_set)),
                     rep("NEG", length(Negative_set))),
                   levels = c("POS", "NEG"))
  
  message("\nFinished creating training data set.")
  
  
  # Select columns that are not informative for ML
  colRem_def = c('chrom','strand','associated_gene', 'associated_transcript', 
             'ref_length','ref_exons', 'ORF_length', 'CDS_length', 'CDS_start', 
             'CDS_end', 'CDS_genomic_start', 'CDS_genomic_end', 'all_canonical',
             'seq_A_downstream_TTS', 'ORF_seq', 'subcategory', 'structural_category',
             'polyA_motif')
  
  # Add RM TP set columns to remove to defined list (if it has been previously created)
  if(exists("colRem_RM")){
    colRem_def <- c(colRem_RM, colRem_def)
  }
  
  # If arg is present, add columns that were provided by --remove_columns
  if(is.null(opt$remove_columns) == FALSE){
    colRem_file <- read.table(opt$remove_columns)
    colRem_file <- unname(unlist(colRem_file))
    
    # Check that all provided names exist
    if(!(all(colRem_file %in% colnames(d)))){
      warning(paste0("\t", "Some columns provided via --remove_columns (-r) do not exist in ",
                     "\n\t", opt$sqanti_classif, " file!"))
    }
    
    # Set both defined and file-input columns as the set to be removed
    colRem <- c(colRem_def,
                colRem_file[which(!(colRem_file %in%  colRem_def))])
    
  }else{
    
    # If arg is not present, only pre-set columns will be removed
    colRem <- colRem_def
  }
 
  # Remove columns from trainingset
  trainingset = trainingset[,-which(colnames(trainingset) %in% colRem)]
  
  
  
  ####  Partition data
  message("\nPartitioning data into training and test sets...")
  message(paste("\n \tProportion of the data to be used for training:", opt$percent_training))
  
  set.seed(123)
  inTraining = caret::createDataPartition(Class, p = opt$percent_training, 
                                          list = FALSE, times = 1)[,1]
  
  training = trainingset[inTraining,]
  testing = trainingset[-inTraining,]
  
  message("\nDescription of the training set:")
  message("\n \tPositive and negative transcript isoforms in training set:")
  print(table(d[rownames(d) %in% rownames(training),]$structural_category))
  message("\n \tPositive and negative transcript isoforms in test set:")
  print(table(d[rownames(d) %in% rownames(testing),]$structural_category))
  
  
  # Train Random Forest classifier with 10 times 10 cross validation
  message("\n-------------------------------------------------")
  
  # check working directory for previously classifier object
  RF_outfiles <- dir(opt$dir)
  
  if("randomforest.RData" %in% RF_outfiles == TRUE){
    
    # if the object exists, it is loaded to save runtime
    
    message(paste0("\nRandom forest classifier already exists in output directory ", 
                  ": loading randomforest.RData object."))
    message("\n\t ***Note: this will skip classifier training.")
    message("\t If you have modified TP and TN sets and wish to train a new model, 
            delete randomforest.RData or provide a different output directory.")
    
    randomforest <- readRDS(paste0(opt$dir, "/randomforest.RData"))
    
  } else {
    
    # if it does not exist, the classifier is trained again
    
    message("\nTraining Random Forest Classifier...")
    message("\n\t***Note: this can take up to several hours.")
    message("\nPre-defined Random Forest parameters (supplied to caret::trainControl()):")
    message("\t - Downsampling in training set (sampling = 'down').")
    message("\t - 10x cross-validation (repeats = 10).\n")
    
    ctrl = caret::trainControl(method = "repeatedcv", repeats = 10,
                      classProbs = TRUE,
                      sampling = 'down', returnData = TRUE,
                      savePredictions = TRUE, returnResamp ='all')
    
    set.seed(1)
    
    randomforest <- caret::train(x = training, 
                                 y = as.factor(Class[inTraining]),
                                 method ='rf',
                                 tuneLength = 15,
                                 metric = "Accuracy",
                                 trControl = ctrl)
    
    saveRDS(randomforest, file = paste(opt$dir, "randomforest.RData",
                                    sep="/"))
    
    message("\nRandom forest training finished.")
    message("\nSaved generated classifier to randomforest.RData file.")
  }
  
  
  #### Apply classifier on test set

  message("\n-------------------------------------------------")
  message("\nRandom forest evaluation: applying classifier to test set...")
  test_pred_prob = predict(randomforest,testing,type = 'prob')
  pred = factor(ifelse(test_pred_prob$POS >= opt$threshold, "POS", "NEG"), 
                levels = c("POS", "NEG"))
  test_result = data.frame(POS = test_pred_prob$POS,
                 NEG = test_pred_prob$NEG,
                 pred, 
                 obs = Class[-inTraining])
  rownames(test_result) <- rownames(test_pred_prob)
  
  # Calculate AUC, sensitivity, specificity
  message("\nTest set evaluation results:")
  message("------------------------------")
  message("\nAUC, Sensitivity and Specificity on test set:")
  print(caret::twoClassSummary(test_result, lev = levels(test_result$obs)))
  
  write.table(data.frame(caret::twoClassSummary(test_result, lev = levels(test_result$obs))),
          file = paste(opt$dir,'/testSet_summary.txt',sep=''),
          quote = F, col.names = F)
  
  message("\nWriting summary to testSet_summary.txt file.")
  
  # Create confusion matrix
  cm = caret::confusionMatrix(data = pred,
                       reference = Class[-inTraining],
                       positive = "POS")
  
  message("\nConfusion matrix:")
  print(cm$table)
  
  message("\nWriting confusion matrix and statistics to output files:")
  
  message("\t testSet_confusionMatrix.txt")
  write.table(data.frame(cm$table),
              file = paste0(opt$dir,"/testSet_confusionMatrix.txt"), 
              quote = F, col.names = T, row.names = F, sep = "\t")
  
  # format stats
  overall <- data.frame(cm$overall)
  colnames(overall) <- "value"
  byClass <- data.frame(cm$byClass)
  colnames(byClass) <- "value"
  cm_stats <- rbind(data.frame(overall),
                    data.frame(byClass))
  
  message("\t testSet_stats.txt")
  write.table(cm_stats, file = paste0(opt$dir,"/testSet_stats.txt"),
              quote = F, col.names = F, sep = "\t")
  
  
  # Variable importance for the prediction
  imp = caret::varImp(randomforest,scale = FALSE)
  imp <- imp$importance
  imp <- data.frame(rownames(imp),imp) 
  imp <- imp[order(-imp$Overall, decreasing = FALSE),]
  imp <- imp[,-1, drop = FALSE]
  
  message("\nGlobal variable importance in Random Forest classifier:")
  print(imp)
  
  write.table(imp, file = paste(opt$dir,
                               "/classifier_variable-importance_table.txt",sep =''), 
              sep = "\t", quote = F, col.names = F)
  
  message("\nVariable importance table saved as classifier_variable-importance_table.txt")
  
  
  ###### ROC curve
  message("\nCalculating and printing test set ROC curves...")
  
  # 1) in function of the probability on the test set: 
  # (not the same proportion of positives and negatives)
  pdf(file = paste0(opt$dir,"/testSet_ROC_curve.pdf"))
  r = pROC::roc(as.numeric(Class[-inTraining]),test_pred_prob$POS,percent = TRUE)
  pROC::auc(r)
  pROC::plot.roc(r, main = "ROC with unbalanced classes")
  text(20,10,paste('AUC =', signif(pROC::auc(r),4)))
  text(20,5,paste('CI 95% = [', signif(pROC::ci(r)[1],4), ',', 
                  signif(pROC::ci(r)[2]),']'))
  
  # 2) same proportion positives and negatives on the test set:
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
  
  r = pROC::roc(as.numeric(Classnewtest),test_pred_prob2$POS,percent = TRUE)
  pROC::auc(r)
  pROC::plot.roc(r, main = "ROC with balanced classes")
  lines(r,col = 'red')
  pROC::ci(r)
  text(20,20,paste('AUC =',signif(pROC::auc(r),4)),col = 'red')
  text(20,15,paste('CI 95% = [',signif(pROC::ci(r)[1],4),',', 
                   signif(pROC::ci(r)[2]),']'),col = 'red')
  dev.off()
  
  message("\nROC curves saved to testSet_ROC_curve.pdf file. Includes:")
  message("\t - ROC curve with unbalanced classes.")
  message("\t - ROC curve with balanced classes.")
  
  
  #################################
  ## Applying classifier to data ##
  #################################
  message("\n------------------------------------------------")
  message("\nApplying Random Forest classifier to input dataset...")
  
  isoform.predict = predict(randomforest,dmult[,colnames(training)],type = 'prob')
  colnames(isoform.predict) = gsub("NEG","NEG_MLprob", colnames(isoform.predict))
  colnames(isoform.predict) = gsub("POS","POS_MLprob", colnames(isoform.predict))
  message("\nRandom forest prediction finished successfully!")
  
  ## Adding predictions to classification table
  classified.isoforms = cbind(dmult[rownames(isoform.predict),],isoform.predict)
  
  negatives = classified.isoforms[classified.isoforms$POS_MLprob < opt$threshold,]
  
  if (opt$force_fsm_in) {
    negatives <- negatives[negatives$structural_category != "full-splice_match",]
  }
  
  classified.isoforms$"ML_classifier" <- "Positive"
  classified.isoforms$"ML_classifier"[rownames(classified.isoforms) %in% 
                                        rownames(negatives)] <- "Negative"
  message("\nRandom forest classification results:")
  print(table(classified.isoforms$ML_classifier))
  
  
  ####################################
  ## Add monoexons back to dataset  ##
  ####################################
    
  # select mono-exons
  dme <- d1[setdiff(rownames(d1), rownames(isoform.predict)),]
  
    # verify that the dme object is not empty 
    # this happens when there are no mono-exons in the input transcriptome
    if(nrow(dme) > 0){
      
      # set ML-related columns to NA
      dme$"NEG_MLprob"  <- dme$"POS_MLprob" <- dme$"ML_classifier"  <- NA
        
      # add mono-exons back to the main data.frame
      d1 <- rbind(classified.isoforms, dme[,colnames(classified.isoforms)])
      
    } else{
      d1 <- classified.isoforms
    }
  
  
# END OF if(run_ML == TRUE)
  
} else if(run_ML == FALSE){
  
  ### Output generation when ML filter is skipped
  
  # ML did not run, fill result columns with NA
  d1$"NEG_MLprob" <- d1$"POS_MLprob" <- d1$"ML_classifier" <- NA
  
}


####################################
###### INTRA-PRIMING FILTERING #####
####################################

message("\n-------------------------------------------------")
message("\nApplying intra-priming filter to our dataset.")

# apply to all not FSM transcripts and to all mono-exons from all categories
d1[,"intra_priming"] <- d1$perc_A_downstream_TTS > as.numeric(opt$intrapriming) & 
  # not FSM: always TRUE
  (d1$structural_category != "full-splice_match" |
     # FSM: only TRUE when mono-exon
     (d1$structural_category == "full-splice_match" & d1$exons == 1))


message("\nIntra-priming filtered transcripts:")
print(table(d1$intra_priming))


# reorder d and d1 to have transcripts in the same order
d1 <- d1[rownames(d),]



#############################################
##### GENERATE AND OUTPUT RESULT TABLES #####
#############################################

message("\n -------------------------------------------------")
message("\nWriting filter results to classification file...")

### Output ML results

# select new columns in d1 that contain the ML results
result_cols <- c("POS_MLprob", "NEG_MLprob", "ML_classifier", "intra_priming")


# add isoform ids and result columns to d object (initial classification table)
ids_df <- data.frame(isoform = rownames(d1))    
d_out <- cbind(ids_df, d, d1[,result_cols])

    
# intersect results of ML and intra-priming to create new column with filter results
d_out$filter_result <- ifelse(d_out$intra_priming == FALSE &
                                (d_out$ML_classifier == "Positive" |
                                   is.na(d_out$ML_classifier)), 
                              yes = "Isoform", no = "Artifact")

    # Condition table:
    # intra_priming == TRUE -> Artifact
    # intra_priming == FALSE:
    #     - ML_classifier == NA -> Isoform
    #     - ML_classifier == Positive -> Isoform
    #     - ML_classifier == Negative -> Artifact


# if mono-exons are to be excluded (--e), set filter_result column to Artifact
if(opt$force_multi_exon == TRUE){
  d_out$filter_result[d_out$exons == 1] <- "Artifact"
}


# write new classification table
write.table(d_out, file = paste0(opt$dir, "/", opt$output, 
                           "_MLresult_classification.txt"),
            quote = FALSE, col.names = TRUE, sep ='\t', row.names = FALSE)

message(paste0("\n\tWrote filter results (ML and intra-priming) to new classification table:\n", 
              "\t", opt$output, "_MLresult_classification.txt file."))



### Generate true isoform list
inclusion_list <- data.frame(Isoforms = rownames(d_out[which(d_out$filter_result == "Isoform"),]))

write.table(inclusion_list, file = paste0(opt$dir, "/", opt$output,
                                    "_inclusion-list.txt"),
            quote = FALSE, col.names = FALSE, sep ='\t', row.names = FALSE)

message(paste0("\n\tWrote isoform list (classified as non-artifacts by both ML and intra-priming", 
               "\n\t", "filters) to ", opt$output, "_inclusion-list.txt file"))



### Output ML-intermediate file if requested

# add isoform ids and write
d1_out <- cbind(ids_df, d1)

if(opt$intermediate_files == TRUE){
  write.table(d1, file = paste0(opt$dir, "/intermediate_", opt$output,
                                "_MLinput_table.txt"),
              quote = FALSE, col.names = TRUE, sep = "\t", row.names = TRUE)
}


# final print summarizing results
message("\n-------------------------------------------------")
message("\nSUMMARY OF MACHINE LEARNING + INTRA-PRIMING FILTERS:")
print(table(d_out$filter_result))

    
message("\n-------------------------------------------------")
message("\nSQANTI3 ML filter finished successfully!")
message("\n-------------------------------------------------")


################################################
