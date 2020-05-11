LR.rarefaction <- function (input, ticks = NULL, rep = 5, samples = 1) {
    if (inherits(input, "eSet") == FALSE) 
        stop("Error. You must give an eSet object\n")
    if (!is.null(assayData(input)$exprs)) {
        data = assayData(input)$exprs
    }
    else {
        data = assayData(input)$counts
    }
    genes <- rownames(data)
    nsamples <- ncol(data)
    if (nsamples > 1) {
        data <- cbind(data, all = rowSums(data))
        samples <- c(samples, nsamples + 1)
    }
    simulation <- myticks <- vector("list", length(samples))
    names(simulation) <- names(myticks) <- colnames(data[, samples])
    for (s in 1:length(samples)) {
        print(paste("sample", s))
        expand.data <- base::rep(genes, data[, samples[s]])
        if (is.null(ticks)) {
            sticks <- round(length(expand.data) * seq(0.1, 1, 
                0.1), 0)
            myticks[[s]] <- sticks
        }
        else {
            myticks[[s]] <- sticks <- ticks
        }
        transcripts <- list()
        for (i in 1:length(sticks)) {
            print(paste("tick", i))
            number = list()
            for (j in 1:rep) {
                print(paste("rep", j))
                number[[j]] <- table(sample(expand.data, sticks[i], 
                  replace = FALSE))
            }
            transcripts[[i]] <- number
        }
        simulation[[s]] <- transcripts
    }
    if (!is.null(featureData(input)$Biotype)) {
        biotype <- as.data.frame(featureData(input)$Biotype)
        rownames(biotype) <- rownames(data)
    }
    if (!is.null(featureData(input)$Category)) {
        category <- as.data.frame(featureData(input)$Category)
        rownames(category) <- rownames(data)
    }
    result <- list(simulation = simulation, ticks = myticks, 
        biotype = biotype, category = category)
}
