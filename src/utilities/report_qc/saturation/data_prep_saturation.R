readData <-
  function (data = NULL, factors = NULL, length = NULL, biotype = NULL, chromosome = NULL, category = NULL, gc = NULL) {
    
    if (is.null(data))
      stop("Expression information must be provided to the readData function")
    
    if (is.null(factors))
      stop("Condition information must be provided to the readData funcion")
    
    if (is.null(length) == FALSE && is.vector(length) == FALSE && is.data.frame(length) == FALSE && is.matrix(length) == FALSE)
      stop( "The length info should be a vector or a data.frame/matrix.")
    
    if (is.null(gc) == FALSE && is.vector(gc) == FALSE && is.data.frame(gc) == FALSE && is.matrix(gc) == FALSE)
      stop( "The GC content info should be a vector or a data.frame/matrix.")
    
    if (is.null(chromosome) == FALSE && ncol(chromosome) != 3)
      stop( "The chromosome object should be a matrix or data.frame with 3 columns: chromosome, start position and end position.")
    
    if (is.null(biotype) == FALSE && is.vector(biotype) == FALSE && is.data.frame(biotype) == FALSE && is.matrix(biotype) == FALSE)
      stop( "The biotype info should be a vector or a data.frame/matrix.")
    
    if (is.null(category) == FALSE && is.vector(category) == FALSE && is.data.frame(category) == FALSE && is.matrix(category) == FALSE)
      stop( "The category info should be a vector or a data.frame/matrix.")
    
    countData <- as.matrix( data )
    
    rowNames <- rownames(countData)
    
    if (nrow(factors) == ncol(countData)) {
      rownames(factors) <- colnames(countData)
    } else {
      stop ("Number of rows in factors must be equal to number of columns in data.\n")
    }
    
    
    pheno <- AnnotatedDataFrame(data=as.data.frame(factors))
    
    input <- ExpressionSet(
      assayData = countData,
      phenoData = pheno)
    
    if (!is.null(length))
      input <- addData(data = input, length = length)
    
    if (!is.null(gc))
      input <- addData(data = input, gc = gc)
    
    if (!is.null(biotype))
      input <- addData(data = input, biotype = biotype)
    
    if (!is.null(chromosome))
      input <- addData(data = input, chromosome = chromosome)
    
    if (!is.null(category))
      input <- addData(data = input, category = category)
    
    input    
    
  }


######################################################
######################################################
######################################################


addData <- function(data, length = NULL, biotype = NULL, chromosome = NULL, factors = NULL, category = NULL, gc = NULL) {
  
  if (inherits(data,"eSet") == FALSE)
    stop("Error. You must give an eSet object.")
  
  if (is.null(length) == FALSE && is.vector(length) == FALSE && is.data.frame(length) == FALSE && is.matrix(length) == FALSE)
    stop( "The length info should be a vector or a data.frame/matrix.")  
  
  if (is.null(gc) == FALSE && is.vector(gc) == FALSE && is.data.frame(gc) == FALSE && is.matrix(gc) == FALSE)
    stop( "The GC content info should be a vector or a data.frame/matrix.")
  
  if (is.null(biotype) == FALSE && is.vector(biotype) == FALSE && is.data.frame(biotype) == FALSE && is.matrix(biotype) == FALSE)
    stop( "The biotype info should be a vector or a data.frame/matrix.")
  
  if (is.null(category) == FALSE && is.vector(category) == FALSE && is.data.frame(category) == FALSE && is.matrix(category) == FALSE)
    stop( "The category info should be a vector or a data.frame/matrix.")
  
  if (is.null(chromosome) == FALSE && ncol(chromosome) != 3)
    stop( "The chromosome object should be a matrix or data.frame with 3 columns: chromosome, start position and end position.")
  
  if (!is.null(assayData(data)$exprs))
    rowNames <- rownames(assayData(data)$exprs)
  else
    rowNames <- rownames(assayData(data)$counts)
  
  # If exists length
  if (!is.null(length)) {
    Length <- rep(NA,length(rowNames))
    names(Length) <- rowNames
    if (is.vector(length)) {
      Length[rowNames] <- as.numeric(as.character(length[rowNames]))
    } else if (is.data.frame(length) || is.matrix(length)) {
      if (ncol(length) == 2) {
        # We assume that the feature names are in the first column and the length in the second
        rownames(length) <- length[,1]
        Length[rowNames] <- as.numeric(as.character( length[rowNames,2]  ))
      } else if (ncol(length) == 1) {
        # We assume that the length are in the first column and the feature names in the rownames
        Length[rowNames] <- as.numeric(as.character( length[rowNames,1]  ))
      } else {
        stop( "The length matrix/data.frame contains more columns than expected.")
      }
    }
    
    featureData(data)@data <- cbind(featureData(data)@data, Length)
  }
  
  # If exists gc
  if (!is.null(gc)) {
    GC <- rep(NA,length(rowNames))
    names(GC) <- rowNames
    if (is.vector(gc)) {
      GC[rowNames] <- as.numeric(as.character(gc[rowNames]))
    } else if (is.data.frame(gc) || is.matrix(gc)) {
      if (ncol(gc) == 2) {
        # We assume that the feature names are in the first column and the GC content in the second
        rownames(gc) <- gc[,1]
        GC[rowNames] <- as.numeric(as.character( gc[rowNames,2]  ))
      } else if (ncol(gc) == 1) {
        # We assume that the GC contents are in the first column and the feature names in the rownames
        GC[rowNames] <- as.numeric(as.character( gc[rowNames,1]  ))
      } else {
        stop( "The GC matrix/data.frame contains more columns than expected.")
      }
    }
    
    featureData(data)@data <- cbind(featureData(data)@data, GC)
  }
  
  # If exists biotype
  if (!is.null(biotype)) {
    Biotype <- rep(NA,length(rowNames))
    names(Biotype) <- rowNames
    if (is.vector(biotype)) {
      Biotype[rowNames] <- as.character(biotype[rowNames])
    } else if (is.data.frame(biotype) || is.matrix(biotype)) {
      if (ncol(biotype) == 2) {
        # We assume that the feature names are in the first column and the biotypes in the second
        rownames(biotype) <- biotype[,1]
        Biotype[rowNames] <- as.character( biotype[rowNames,2]  )
      } else if (ncol(biotype) == 1) {
        # We assume that the biotypes are in the first column and the feature names in the rownames
        Biotype[rowNames] <- as.character( biotype[rowNames,1]  )
      } else {
        stop( "The biotype matrix/data.frame contains more columns than expected.")
      }
    }
    
    featureData(data)@data <- cbind(featureData(data)@data, Biotype)
    featureData(data)@data$Biotype <- as.character(featureData(data)@data$Biotype)
  } 
  
  
  # If exists biotype
  if (!is.null(category)) {
    Category <- rep(NA,length(rowNames))
    names(Category) <- rowNames
    if (is.vector(category)) {
      Category[rowNames] <- as.character(category[rowNames])
    } else if (is.data.frame(category) || is.matrix(category)) {
      if (ncol(category) == 2) {
        # We assume that the feature names are in the first column and the biotypes in the second
        rownames(category) <- category[,1]
        Category[rowNames] <- as.character( category[rowNames,2]  )
      } else if (ncol(category) == 1) {
        # We assume that the biotypes are in the first column and the feature names in the rownames
        Category[rowNames] <- as.character( category[rowNames,1]  )
      } else {
        stop( "The category matrix/data.frame contains more columns than expected.")
      }
    }
    
    featureData(data)@data <- cbind(featureData(data)@data, Category)
    featureData(data)@data$Category <- as.character(featureData(data)@data$Category)
  } 
  
  # If exists chromosome
  if (!is.null(chromosome)) {
    Chromosome <- GeneStart <- GeneEnd <- rep(NA,length(rowNames))
    names(Chromosome) <- names(GeneStart) <- names(GeneEnd) <- rowNames
    
    Chromosome[rowNames] <- as.character(chromosome[rowNames,1])
    GeneStart[rowNames] <- as.numeric(as.character(chromosome[rowNames,2]))
    GeneEnd[rowNames] <- as.numeric(as.character(chromosome[rowNames,3]))
    
    featureData(data)@data <- cbind(featureData(data)@data, Chromosome, GeneStart, GeneEnd)
  }
  
  # If exists new factors
  if (!is.null(factors))
    phenoData(data)@data <- cbind(phenoData(data)@data, factors)
  
  data
  
}
