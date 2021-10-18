#######################################
#                                     #
#  SQANTI3 output comparison report   #
#             generation              #
#                                     #
#######################################


# Author: Jorge Martinez Tomas & Alejandro Paniagua
# Last modified: 21/09/2021 by Jorge Martinez


#######################################
#                                     #
#      PACKAGES AND LIBRARIES         #
#                                     #
#######################################

suppressMessages(library(DT))
suppressMessages(library(gridExtra))
suppressMessages(library(knitr))
suppressMessages(library(optparse))
suppressMessages(library(rmarkdown))
suppressMessages(library(tidyverse))
suppressMessages(library(UpSetR))
suppressMessages(library(VennDiagram))


#######################################
#                                     #
#             FUNCTIONS               #
#                                     #
#######################################

# -------------------- Read data

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# -------------------- Tags and basic comparison P/A

count_FL <- function(count_f, class_f){
  # Count FL reads from TSV file
  count_f <- data.table::data.table(count_f)
  count_f.compact <- count_f[, list(
    value=length(read_id)
  ), by="transcript_id"]
  count_f.compact <- as.data.frame(count_f.compact)
  names(count_f.compact)[1] <- "isoform" 
  
  count_f.join <- list(class_f, count_f.compact) %>% 
    purrr::reduce(left_join, by="isoform")
  class_f$FL <- count_f.join$value
  return(class_f)
}

isoformTags <- function(junctions_file) {
  # Create unique junction chains (UJC)
  # Builds isoform tags (UJC) with the following structure: Chr_strand_start_end
  df <- junctions_file[, c("isoform", "chrom", "strand")] # df with isoforms in *junctions.txt
  dt <- data.table::data.table(df)
  dt <- dt[,coord:=paste0(junctions_file$genomic_start_coord, "_", junctions_file$genomic_end_coord)]
  dt <-
    dt[, list(tagcoord = paste0(coord, collapse = "_")),
       by = c("isoform", "chrom", "strand")]
  df <- as.data.frame(dt)
  df$tags <- paste(df$chrom, df$strand, df$tagcoord, sep = "_")
  #df <- df[order(df$isoform),]
  #tag <- paste(df$chrom, df$strand, df$tagcoord, sep = "_")
  return(df[,c("isoform","tags")])
}


addSC <- function(class_file){
  # Add the structural category to the tag (UJC)
  str_cat <- c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog")
  shortSC <- c("FSM", "ISM", "NIC", "NNC")
  
  class_file$SC <- class_file$structural_category
  for (i in 1:length(str_cat)){
    cond <- which(class_file$SC == str_cat[i])
    class_file$SC[cond] <- shortSC[i]
  }
  class_file$tags <- paste0(class_file$SC, "_", class_file$tags)
  class_file$SC <- NULL
  return(class_file)
}


filter_monoexon <- function(class_file){
  # Delete monoexons
  filtered_classification <- class_file[class_file$exons > 1, ]
  filtered_classification[order(filtered_classification$isoform),]
}


filter_monoexon_sirv <- function(class_file){
  # Delete monoexons and transcripts from rare chromosomes
  valid_chrom <- c(paste0("chr", 1:22), paste0("chr", c("X","Y")), paste0("SIRV", 1:7))
  
  filtered_classification <- class_file[class_file$exons > 1 &
                                          class_file$chrom %in% valid_chrom, ]
  filtered_classification[order(filtered_classification$isoform),]
}


swapcoord <- function(dfclass){
  # Swap TSS and TTS coordinates depending on the strand
  tss <- dfclass$TTS_genomic_coord[dfclass$strand == "-"]
  tts <- dfclass$TSS_genomic_coord[dfclass$strand == "-"]
  dfclass$TSS_genomic_coord[dfclass$strand == "-"] <- tss
  dfclass$TTS_genomic_coord[dfclass$strand == "-"] <- tts
  return(dfclass)
}


reverseswap <- function(class_file) {
  # Reverse the swap from swapcoord()
  class_file <- swapcoord(class_file)
  sdtts <- class_file$sdTSS[class_file$strand == "-"]
  sdtss <- class_file$sdTTS[class_file$strand == "-"]
  
  class_file$sdTSS[class_file$strand == "-"] <- sdtss
  class_file$sdTTS[class_file$strand == "-"] <- sdtts
  return(class_file)
}


uniquetag <- function(class_file) {
  # Aggregate by UJC and calculate many metrics when coords ar given
  dt <- data.table::data.table(class_file)
  dt.out <-
    dt[, list(
      associated_gene=list(unique(associated_gene)),
      FL=as.numeric(sum(FL)),
      l_FL=list(FL),
      exons=unique(exons),
      TSS_genomic_coord=list(TSS_genomic_coord),
      TTS_genomic_coord=list(TTS_genomic_coord),
      length=list(length),
      medianTSS = round(median(TSS_genomic_coord)),
      medianTTS = round(median(TTS_genomic_coord)),
      sdTSS = sd(TSS_genomic_coord),
      sdTTS = sd(TTS_genomic_coord),
      isoform = list(isoform)
    ), by = c("tags", "structural_category")]
  dt.out <- as.data.frame(dt.out)
  return(dt.out[order(dt.out$tags, dt.out$structural_category),])
}


uniquetag_simple <- function(class_file) {
  # Aggregate by UJC and calculate some basic metrics
  dt <- data.table::data.table(class_file)
  dt.out <-
    dt[, list(
      associated_gene=list(unique(associated_gene)),
      FL=as.numeric(sum(FL)),
      l_FL=list(FL),
      exons=unique(exons),
      length=list(length),
      isoform = list(isoform)
    ), by = c("tags", "structural_category")]
  dt.out <- as.data.frame(dt.out)
  return(dt.out[order(dt.out$tags, dt.out$structural_category),])
}


multipleComparison <- function(l_class){
  # Presence and ausence of UJC in the samples
  a <- c(rbind(names(l_class), paste0(names(l_class), "SC")))
  
  l_class %>%
    purrr::map(~ data.frame(col = .$tags, .$tags,.$structural_category, stringsAsFactors = FALSE)) %>%
    purrr::reduce(full_join, by = "col") %>%
    select(-col) %>%
    setNames(a)
}


# -------------------- Isoform (UJC) analysis

SD_TSS_TTS <- function(l_class){
  # Calculate UJC SD
  a <- c("tags", rbind(paste0(names(l_class), "TSS"), paste0(names(l_class), "TTS")))
  TSS_TTS_params <- list()
  for (i in 1:length(l_class)){
    TSS_TTS_params[[i]] <- l_class[[i]][,c("tags", "TSS_genomic_coord", "TTS_genomic_coord")]
  }
  TSS_TTS_params <- TSS_TTS_params %>% 
    purrr::reduce(full_join, by="tags") %>% 
    setNames(a)
  
  a <- paste0(names(l_class), "TSS")
  b <- paste0(names(l_class), "TTS")
  allTSS <- TSS_TTS_params[, a]
  allTSS[allTSS == "NULL"] <- NA
  allTTS <- TSS_TTS_params[, b]
  allTTS[allTTS == "NULL"] <- NA
  
  TSS_TTS_df <- data.frame(tags=TSS_TTS_params$tags)
  
  # max and min value
  sapplycolumns <- function(data, func){
    tmp <- list()
    for (i in 1:ncol(data)){
      tmp[[names(data)[i]]] <- sapply(data[,i], func)
    }
    for (i in 1:length(tmp)){
      tmp[[i]][tmp[[i]]=="NULL"] <- NA
    }
    return(as.data.frame(tmp))
  }
  
  minNA <- function(x) ifelse(length(x) > 1, min(x), NA)
  
  
  maxTSS <- sapplycolumns(allTSS, max)
  maxTSS[maxTSS=="-Inf"] <- NA
  
  maxTTS <- sapplycolumns(allTTS, max)
  maxTTS[maxTTS=="-Inf"] <- NA
  
  minTSS <- sapplycolumns(allTSS, minNA)
  minTSS[minTSS=="Inf"] <- NA
  
  minTTS <- sapplycolumns(allTTS, minNA)
  minTTS[minTTS=="Inf"] <- NA
  
  minmaxTSS <- cbind(minTSS, maxTSS)
  minmaxTTS <- cbind(minTTS, maxTTS)
  
  TSS_TTS_df$minmax.SD.TSS <- apply(minmaxTSS, 1, function(x) sd(unlist(x), na.rm = TRUE))
  
  TSS_TTS_df$minmax.SD.TTS <- apply(minmaxTTS, 1, function(x) sd(unlist(x),na.rm = TRUE))
  
  # Median value
  
  medianTSS <- sapplycolumns(allTSS, median)
  medianTTS <- sapplycolumns(allTTS, median)
  
  TSS_TTS_df$median.SD.TSS <- apply(medianTSS, 1, function(x) sd(unlist(x),na.rm = TRUE))
  
  TSS_TTS_df$median.SD.TTS <- apply(medianTTS, 1, function(x) sd(unlist(x),na.rm = TRUE))
  
  # Max SD
  
  TSS_TTS_df$SD.TSS <- apply(TSS_TTS_df[,c("minmax.SD.TSS", "median.SD.TSS")],1,max)
  TSS_TTS_df$SD.TTS <- apply(TSS_TTS_df[,c("minmax.SD.TTS", "median.SD.TTS")],1,max)
  
  return(TSS_TTS_df[, c("tags", "SD.TSS", "SD.TTS")])
}

iso_analysis <- function(l_class){
  
  class_bind <- bind_rows(l_class)
  class_bind <- data.table::data.table(class_bind)
  
  class_compact <- class_bind[, list(
    exons=unique(exons),
    length=as.numeric(median(unlist(length))),
    FL_cpm=as.numeric(median(FL))
  ), by="tags"]
  
  class_compact <- as.data.frame(class_compact)
  tot_reads <- sum(class_compact$FL_cpm)
  
  class_compact$FL_cpm <- (class_compact$FL_cpm * 10^6)/tot_reads 
  
  if (TSS_TTS_coord == TRUE) {
    TSS_TTS_df <- SD_TSS_TTS(l_class)
    df_iso <- list(TSS_TTS_df, class_compact) %>% 
      purrr::reduce(full_join, by="tags")
  } else {df_iso <- data.frame(
    tags = class_compact$tags,
    SD.TSS=NA,
    SD.TTS=NA,
    exons=class_compact$exons,
    length=class_compact$length,
    iso_exp=class_compact$iso_exp
  )}
  
  return(df_iso)
  
}


# -------------------- Gene analysis

del_novel_genes <- function(df){
  df.out <- df[!grepl("novelGene|SIRV", df$associated_gene), ]
  df.out$associated_gene <- as.character(df.out$associated_gene)
  return(df.out)
}

tags_per_gene <- function(res_class){
  # Delete novel genes and SIRV
  ref_genes <- del_novel_genes(res_class)
  ref_genes <- data.table::data.table(ref_genes)
  
  ref_genes.out <- ref_genes[, list(
    N_UJC=length(unique(tags)),
    tags=list(tags)
  ), by="associated_gene"]
  
  return(as.data.frame(ref_genes.out))
}

jaccard <- function(x, comb){
  x <- x[2:length(x)]
  val <- c()
  for (i in 1:ncol(comb)){
    inter <- length(intersect(x[[comb[1,i]]], x[[comb[2,i]]]))
    uni <- length(union(x[[comb[1,i]]], x[[comb[2,i]]]))
    jac <- inter/uni
    val <- c(val, jac)
  }
  return(median(val))
}

get_jaccard <- function(l_df, l_class){
  a <- c("associated_gene", paste0("tags_", names(l_class)))
  l <- list()
  for (i in 2:length(l_df)){
    l[[i-1]] <- l_df[[i]][,c("associated_gene","tags")]
  }
  df <- l %>% 
    purrr::reduce(full_join, by="associated_gene") %>% 
    setNames(a)
  
  n <- colnames(df[,2:ncol(df)])
  comb <- combn(n, 2)
  
  jc <- apply(df, 1, function(x) jaccard(x, comb))
  df <- data.frame(associated_gene=df$associated_gene, jaccard=jc)
  return(df)
}

gene_expr <- function(res_class){
  ref_genes <- del_novel_genes(res_class)
  ref_genes <- data.table::data.table(ref_genes)
  
  ref_genes.UJC <- ref_genes[, list(
    associated_gene=unique(associated_gene),
    FL=as.numeric(median(FL))
  ), by="tags"]
  
  ref_genes.out <- ref_genes.UJC[, list(
    FL_cpm=sum(FL)
  ), by="associated_gene"]
  
  ref_genes.out <- as.data.frame(ref_genes.out)
  
  tot_exp <- sum(ref_genes.out$FL_cpm)
  ref_genes.out$FL_cpm <- (ref_genes.out$FL_cpm * 10^6) / tot_exp
  
  return(ref_genes.out)
}

gene_length <- function(df_gtf, df_gene){
  cond <- which(df_gtf$feature == "exon" & df_gtf$associated_gene %in% df_gene$associated_gene)
  df_gtf <- df_gtf[cond,]
  df_gtf$diff <- df_gtf$end - df_gtf$start
  df_gtf <- data.table::data.table(df_gtf)
  df_gtf.out <- df_gtf[, list(
    length=sum(diff)
  ), by="associated_gene"]
  return(as.data.frame(df_gtf.out))
}

gene_analysis <- function(l_class){
  l_res_class <- list()
  for (i in 1:length(l_class)){
    l_res_class[[i]] <- del_novel_genes(l_class[[i]])
  }
  class_bind <- bind_rows(l_res_class)
  
  l_gene_df <- list()
  l_gene_df[[1]] <- tags_per_gene(class_bind)
  for (i in 1:length(l_class)){
    l_gene_df[[i+1]] <- tags_per_gene(l_class[[i]])
  }
  
  a <- c("associated_gene", "N_UJC", paste0("N_UJC_", names(l_class)))
  df_gene <- l_gene_df %>% 
    purrr::map(~ data.frame(col = .$associated_gene, .$N_UJC, stringsAsFactors = FALSE)) %>% 
    purrr::reduce(full_join, by="col") %>% 
    setNames(a)
  
  df_jac <- get_jaccard(l_gene_df, l_class)
  
  df_exp <- gene_expr(class_bind)
  if (class(ref_gtf) == "data.frame"){
    df_len <- gene_length(ref_gtf, df_gene)
    l_df_metrics <- list(df_gene, df_jac, df_exp, df_len)
  } else{
    l_df_metrics <- list(df_gene, df_jac, df_exp)
  }
  
  
  df_gene <- l_df_metrics %>% 
    purrr::reduce(full_join, by="associated_gene")
  
  return(df_gene)
}


# -------------------- Final comparison function

compareTranscriptomes <- function(l_iso){
  
  # Add FL count
  if (count.tsv == TRUE){
    print("Counting FL reads...")
    for (i in 1:length(l_iso)){
      l_iso[[i]][[1]] <- count_FL(count.res[[i]] ,l_iso[[i]][[1]])
    }
  }
  
  # Filter monoexon, add tag column and aggregate
  print("Generating unique junction chains (UJC)...")
  n <- names(l_iso)
  l_class <- list()
  for ( i in 1:length(l_iso)) {
    class.filtered <- filter_monoexon(l_iso[[i]][[1]]) # delete monoexons
    
    tags <- isoformTags(l_iso[[i]][[2]]) # build tags
    if (length(tags) != nrow(class.filtered)){
      class.filtered <- filter_monoexon_sirv(l_iso[[i]][[1]])
    }
    #class.filtered[,"tags"] <- tags # add tags
    class.filtered <- list(class.filtered, tags) %>% 
      purrr:::reduce(left_join, by="isoform")
    
    if (TSS_TTS_coord){
      #class.swap <- swapcoord(class.filtered) # swap - strand
      #class.uniquetags <- as.data.frame(uniquetag(class.swap)) # group by tag
      class.out <- as.data.frame(uniquetag(class.filtered)) 
    } else{
      class.out <- uniquetag_simple(class.filtered) # group by tag
    }
    
    class.out <- addSC(class.out)
    l_class[[i]] <- class.out # add to list
  }
  
  print("Performing presence/ausence comparison...")
  names(l_class) <- n # add names
  comptags <- multipleComparison(l_class)
  comptags.SC <- cbind(structural_category =
                         do.call(dplyr::coalesce, comptags[,paste0(n,"SC")]),
                       comptags[,n]
  )
  comptags.out <- cbind( TAGS =
                           do.call(dplyr::coalesce, comptags[,n]),
                         comptags.SC
  )
  
  comptags.PA <- comptags.out
  comptags.PA[,3:ncol(comptags.PA)][!is.na(comptags.PA[,3:ncol(comptags.PA)])] <- 1
  comptags.PA[,3:ncol(comptags.PA)][is.na(comptags.PA[,3:ncol(comptags.PA)])] <- 0
  
  print("Performing UJC analysis...")
  iso.metrics <- iso_analysis(l_class)
  
  print("Analysing associated genes...")
  gene.metrics <- gene_analysis(l_class)
  
  output <-
    list(
      classifications = l_class,
      comparison = comptags.out,
      comparisonPA = comptags.PA,
      iso_metrics = iso.metrics,
      gene_metrics = gene.metrics
    )
  
  return(output)
}


# -------------------- Functions for the report generation

iso2url <- function(id){
  # Isoform coords to Genome Browser
  if (id != "NA") {
    id.split <- str_split(id,"_")
    chr <- id.split[[1]][2]
    start <- id.split[[1]][4]
    end <- id.split[[1]][length(id.split[[1]])]
    name <- substr(id, 1, 30)
    url <- paste0(
      "<a href='https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=",
      chr, "%3A", start, "%2D", end, "&hgsid=1143169919_jAraPbUWtMCdAHfgTHk4sDHQHW7R",
      "'>", name, "...</a>")
    return(url)
  } else { return(NA) }
}



#######################################
#                                     #
#                MAIN                 #
#                                     #
#######################################

# -------------------- Argument parser

option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = NULL,
              help="directory with input files (classification and junction files)",
              metavar = "DIRIN"),
  make_option(c("-o", "--outdir"), type = "character", default = ".",
              help="Output directory for the report and CSV file [default= %default]",
              metavar = "DIROUT"),
  make_option(c("-n", "--name"), type = "character", default = "comparison_output",
              help="Output name for the HTML report and CSV file (without extension) [default= %default]",
              metavar = "OUTNAME"),
  make_option(c("--lrgasp"), action="store_true",type = "character", default = FALSE,
              help="Use lrgasp metrics",
              metavar = "LRGASP")
)

opt_parser = OptionParser(
  usage = "usage: %prog [-i DIRIN] [-o DIROUT] [-n OUTNAME] [--lrgasp]",
  option_list=option_list
)

opt = parse_args(opt_parser)

directory <- opt$dir
output_directory <- opt$outdir
output_name <- opt$name
lrgasp <- opt$lrgasp

if (is.null(directory)) {
  stop("\n\nAt least one argument must be supplied.\nThe -d argument is required (directory containing input files)")
}


# -------------------- Load data

print("Reading input data...")
if (dir.exists(directory)){
  dir_in <- directory
} else {
  dir_in <- paste(getwd(), directory, sep="/")
  if (!dir.exists(dir_in)){
    stop(paste0("\n\nCould not find the input directory (", directory, ").\nPlease enter a valid path"))
  }
}

class_in <-
  list.files(dir_in,
             pattern = "*_classification.txt",
             all.files = FALSE,
             full.names = TRUE)
junct_in <-
  list.files(dir_in,
             pattern = "*_junctions.txt",
             all.files = FALSE,
             full.names = TRUE)

if (length(class_in) != length(junct_in)){
  stop("ERROR: There is a different number of classification and junction files in the directory")
} else if (length(class_in) == 0){
  stop(paste0("ERROR: No classification and junction files were found in the directory: ", dir_in))
}

f_in <- list()
for (i in 1:length(class_in)) {
  f <- class_in[[i]]
  start <- stringr::str_locate(f, dir_in)[[2]]
  end <- stringr::str_locate(f, "_classification.txt")[[1]]
  idx <- substring(f, (start+2), (end-1))
  classification <- read.table(class_in[[i]], header = T, sep = "\t")
  junctions <- read.table(junct_in[[i]], header = T, sep = "\t")
  f_in[[idx]] <- list(classification, junctions)
}

# ----- TSV Count input

count.files <- 
  list.files(dir_in,
             pattern = "*.tsv",
             all.files = FALSE,
             full.names = TRUE)

count.tsv <- TRUE
if (length(class_in) != length(count.files)){
  print("ERROR: Issue loading count files (.tsv)")
  print("Different number of count files (.tsv) than samples")
  count.tsv <- FALSE
} else {
  count.res <- list()
  for (i in 1:length(count.files)){
    f <- count.files[i]
    count.res[[names(f_in)[[i]]]] <- read.table(f, header=TRUE, sep="\t")
  }
}


# ----- GTF input

gtf_name <- list.files(dir_in,
                       pattern = "*reference.gtf",
                       all.files = FALSE,
                       full.names = TRUE)
if (length(gtf_name) == 1){
  ref_gtf <- try({
    f1 <- "first_half_gtf_sqanti_comparator.txt"
    f2 <- "second_half_gtf_sqanti_comparator.txt"
    system(paste0('cut -f3-5 ', gtf_name[[1]], ' > ', dir_in, '/', f1))
    system(paste0('cut -f9 ', gtf_name[[1]], ' | cut -d ";" -f1 |grep -v "#" |cut -d " " -f2 > ', dir_in,'/', f2))
    gtf1 <- read.table(paste0(dir_in, '/', f1), header=FALSE, sep="\t")
    gtf2 <- read.table(paste0(dir_in,'/', f2), header=FALSE, sep="\t", quote = '"')
    full_gtf <- cbind(gtf1, gtf2)
    colnames(full_gtf) <- c("feature", "start", "end", "associated_gene")
    system(paste0('rm ', dir_in, '/', f1))
    system(paste0('rm ', dir_in, '/', f2))
    full_gtf
  }, silent = TRUE)
  if (class(ref_gtf) == "try-error"){
    print("ERROR: Issue loading reference GTF file")
  }
} else{
  print("ERROR: Issue loading reference GTF file")
}


# ----- LRGASP input

if (lrgasp == TRUE){
  lrgasp.files <- 
    list.files(dir_in,
               pattern = "*_results.RData",
               all.files = FALSE,
               full.names = TRUE)
  if (length(class_in) != length(lrgasp.files)){
    print("ERROR: Issue loading LRGASP files")
    print("Different number of LRGASP files than samples")
    lrgasp <- FALSE
  } else {
    lrgasp.res <- list()
    for (i in 1:length(lrgasp.files)){
      f <- lrgasp.files[i]
      lrgasp.res[[names(f_in)[[i]]]] <- loadRData(f)
    }
  }
}


# -------------------- Check for TSS and TTS genomic coords

TSS_TTS_coord <- TRUE
for ( i in 1:length(f_in)){
  if (!("TSS_genomic_coord" %in% colnames(f_in[[i]][[1]]))){
    TSS_TTS_coord <- FALSE
  }
}

# --------------------  P/A comparison and Iso & Gene analysis

res <- try({
  compareTranscriptomes(f_in)
}, silent = TRUE)

if (class(res) == "try-error"){
  print("ERROR: An error has ocurred during the comparison. A partial report will be generated. Here's the original error message:")
  print(geterrmessage())
}


# -------------------- Comparison with res classification files
print("Comparing comparisson results between samples...")

# ----- Define max number of samples in plots

if (length(f_in) <= 6){
  limit <- length(f_in)
} else {limit <- 6}


# ----- Vector of structural categories

#str_cat <- unique(f_in[[1]][[1]]$structural_category)
#str_cat <- c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "antisense", "fusion", "genic", "intergenic")
str_cat <- c("FSM", "ISM", "NIC", "NNC", "Antisense", "Fusion", "Genic-Genomic", "Genic-Intron", "Intergenic")


# ----- Generates dataframe with a summary of the SQANTI3 classification files

df_summary.1 <- data.frame(ID = names(f_in))

# Count total isoform from SQANTI3
n <- c()
for (i in f_in) {
  n <- c(n, nrow(i[[1]]))
}
df_summary.1$total <- n

# Count isoforms for each structural category
for (i in str_cat) {
  n <- c()
  for (j in f_in) {
    k <- j[[1]]
    n <- c(n, nrow(k[k$structural_category == i,]))
  }
  df_summary.1[,i] <- n
}


# ----- Generates dataframe with a summary of the unique tag comparison

df_summary.2 <- data.frame(ID = names(res[[2]][3:ncol(res[[2]])]))

# Count unique tags
n <- c()
for (i in res[[1]]) {
  n <- c(n, nrow(i))
}
df_summary.2$uniq_id <- n

# Count unique tags for each category
for (i in str_cat) {
  n <- c()
  for (j in res[[1]]) {
    n <- c(n, nrow(j[j$structural_category == i,]))
  }
  df_summary.2[,i] <- n
}


# ----- Add GenomeBrowser URL to the P/A table

df.PA <- res[[3]]
df.PA$TAGS <- lapply(df.PA$TAGS, iso2url)


# ----- Counts per gene and exon structure

countpergene <- c()
exonstructure <- c()
for (i in 1:limit){
  data <- f_in[[i]]
  data.class <- data[[1]]
  
  data.class$novelGene <- "Annotated Genes"
  data.class[grep("novelGene", data.class$associated_gene), "novelGene"] <- "Novel Genes"
  data.class$novelGene = factor(data.class$novelGene,
                                levels = c("Novel Genes","Annotated Genes"),
                                ordered=TRUE)
  
  isoPerGene = aggregate(data.class$isoform,
                         by = list("associatedGene" = data.class$associated_gene,
                                   "novelGene" = data.class$novelGene,
                                   "FSM_class" = data.class$FSM_class),
                         length)
  
  data.class[which(data.class$exons>1), "exonCat"] <- "Multi-Exon"
  data.class[which(data.class$exons==1), "exonCat"] <- "Mono-Exon"
  data.class$exonCat = factor(data.class$exonCat,
                              levels = c("Multi-Exon","Mono-Exon"),
                              ordered=TRUE)
  
  canonical.labels=c("Canonical", "Non-canonical")
  data.class$all_canonical = factor(data.class$all_canonical,
                                    labels=canonical.labels,
                                    levels = c("canonical","non_canonical"),
                                    ordered=TRUE)
  
  countpergene <- c(
    countpergene,
    sum(isoPerGene$x == 1),
    sum(isoPerGene$x == 2 | isoPerGene$x == 3),
    sum(isoPerGene$x == 4 | isoPerGene$x == 5),
    sum(isoPerGene$x >= 6)
  )
  
  exonstructure <- c(
    exonstructure,
    sum(data.class$novelGene == "Novel Genes" & data.class$exonCat == "Mono-Exon"),
    sum(data.class$novelGene == "Novel Genes" & data.class$exonCat == "Multi-Exon"),
    sum(data.class$novelGene == "Annotated Genes" & data.class$exonCat == "Mono-Exon"),
    sum(data.class$novelGene == "Annotated Genes" & data.class$exonCat == "Multi-Exon")
    
  )
}

sample <- c(rep(names(f_in[1:limit]), each=4))
number <- rep(c("1","2-3","4-5", ">=6"), times=limit)
isoPerGene <- data.frame(sample, number, countpergene)

category <- rep(c("Novel-Mono", "Novel-Multi", "Annotated-Mono", "Annotated-Multi"), times=limit)
exonstructure <- data.frame(sample, category, exonstructure)


# ----- Summary dataframe pivoted

df_SC <- df_summary.1[1:limit,]
df_SC$total <- NULL
df_SC <- df_SC %>% 
  pivot_longer(!"ID", "SC")

# ----- Distance to TSS, TTS and CAGE peak

dist.list <- list()
dist.msr <- c("diff_to_TSS", "diff_to_TTS", "dist_to_cage_peak")
dist.SC <- c("FSM", "ISM")
contador <- 1
for (i in dist.msr){
  for (j in dist.SC){
    sample <- c()
    dist <- c()
    for (k in 1:limit){
      data.class <- f_in[[k]][[1]]
      cond <- which(data.class$structural_category == j)
      x <- data.class[cond, i]
      sample <- c(
        sample,
        rep(names(f_in)[k], times=length(x))
      )
      dist <- c(
        dist,
        x
      )
    }
    dist.list[[contador]] <- data.frame(sample, dist)
    names(dist.list)[length(dist.list)] <- paste(i,j,sep="_")
    contador <- contador + 1
  }
}

# ----- RT-switching

FSM <- c()
NIC <- c()
NNC <- c()
for (i in 1:limit){
  data.class <- f_in[[i]][[1]]
  df <- group_by(data.class, structural_category, RTS_stage) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
  FSM.match <- df$count[which(df$structural_category == "FSM")]
  FSM <- c(FSM, ((FSM.match[2]/(FSM.match[1]+FSM.match[2]))*100))
  
  NIC.match <- df$count[which(df$structural_category == "NIC")]
  NIC <- c(NIC, ((NIC.match[2]/(NIC.match[1]+NIC.match[2]))*100))
  
  NNC.match <- df$count[which(df$structural_category == "NNC")]
  NNC <- c(NNC, ((NNC.match[2]/(NNC.match[1]+NNC.match[2]))*100))
}

sample <- names(f_in)[1:limit]
FSM.RT <- data.frame(sample, FSM)
NIC.RT <- data.frame(sample, NIC)
NNC.RT <- data.frame(sample, NNC)

# ----- List of unique tags for each sample

l <- list()
for (i in 3:ncol(res[[2]])){
  l[[i-2]] <- na.omit(res[[2]][,i])
}
names(l) <- colnames(res[[2]])[3:ncol(res[[2]])]


# -------------------- LRGASP metrics

if (lrgasp == TRUE){
  sirv.metrics <- list()
  for (i in 1:length(lrgasp.res)){
    sirv.metrics[[i]] <- lrgasp.res[[names(lrgasp.res)[i]]]["SIRV"]
  }
  sirv.metrics <- bind_cols(sirv.metrics)
  colnames(sirv.metrics) <- names(lrgasp.res)
  # FSM
  fsm.metrics <- list()
  for (i in 1:length(lrgasp.res)){
    fsm.metrics[[i]] <- lrgasp.res[[names(lrgasp.res)[i]]]["FSM"]
  }
  fsm.metrics <- bind_cols(fsm.metrics)
  colnames(fsm.metrics) <- names(lrgasp.res)

  # ISM
  ism.metrics <- list()
  for (i in 1:length(lrgasp.res)){
    ism.metrics[[i]] <- lrgasp.res[[names(lrgasp.res)[i]]]["ISM"]
  }
  ism.metrics <- bind_cols(ism.metrics)
  colnames(ism.metrics) <- names(lrgasp.res)

  # NIC
  nic.metrics <- list()
  for (i in 1:length(lrgasp.res)){
    nic.metrics[[i]] <- lrgasp.res[[names(lrgasp.res)[i]]]["NIC"]
  }
  nic.metrics <- bind_cols(nic.metrics)
  colnames(nic.metrics) <- names(lrgasp.res)

  # NNC
  nnc.metrics <- list()
  for (i in 1:length(lrgasp.res)){
    nnc.metrics[[i]] <- lrgasp.res[[names(lrgasp.res)[i]]]["NNC"]
  }
  nnc.metrics <- bind_cols(nnc.metrics)
  colnames(nnc.metrics) <- names(lrgasp.res)
}


#######################################
#                                     #
#     TABLE AND PLOT GENERATION       #
#                                     #
#######################################

# -------------------- 
# -------------------- 
# TABLE INDEX
# t1: summary table
# t2: presence/ausence table
# t3: SIRV metrics

# -------------------- 
# -------------------- 
# PLOT INDEX
# p1: gene characterization
#   p1.1: isoforms per gene
#   p1.2: exon structure
# p2: structural category distribution
# p3: splice junction distribution for each SC
# p4: distance to TSS
# p5: distance to TTS
# p6: distance to CAGE peak
# p7: bad quality features
# p8: good quality features
# p9: Venn diagrams 
# p10: UpSet plot
# p11: Venn diagrams for SC
# p12: UpSet plot for SC
# p13: TSS standard deviation
# p14: TTS standard deviation
# p16: iso analysis
# p17: gene analysis


print("Generating plots for the report...")

# -------------------- Global plot parameters
# COPY-PASTE FROM SQANTI3 REPORT

#myPalette = c("#6BAED6","#FC8D59","#78C679","#EE6A50","#969696","#66C2A4", "goldenrod1", "darksalmon", "#41B6C4","tomato3", "#FE9929")
myPalette = c("#6BAED6","#FC8D59","#78C679","#EE6A50","#969696","#66C2A4", "goldenrod1", "darksalmon", "#41B6C4","tomato3", "#FE9929")

mytheme <- theme_classic(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=13),
        axis.text.x  = element_text(size=12),
        axis.title.y = element_text(size=13),
        axis.text.y  = element_text(vjust=0.5, size=12) ) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(size=12), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=15, hjust = 0.5)) +
  theme(plot.margin = unit(c(2.5,1,1,1), "cm"))

# -------------------- 
# TABLE 1: summary table

t1.1 <- DT::datatable(df_summary.1,
              extensions = "Buttons",
              options = list(
                order = list(list(5, "asc"),list(1, "desc")),
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'pdf', 'print')
              ),
              escape = FALSE,
              caption = "Table 1. Summary from SQANTI3 output comparison")

t1.2 <- DT::datatable(df_summary.2,
                      extensions = "Buttons",
                      options = list(
                        order = list(list(5, "asc"),list(1, "desc")),
                        dom = 'Bfrtip',
                        buttons = c('copy', 'csv', 'pdf', 'print')
                      ),
                      escape = FALSE,
                      caption = "Table 2. Summary from isoform id (tag) comparison")

# TABLE 2: presence/ausence isoforms

t2 <- DT::datatable(df.PA,
              escape = FALSE,
              options = list(
                pageLength = 10,
                autoWidth = TRUE,
                columnDefs = list(list(width = '10px', targets = "_all"))
              ),
              rownames = FALSE,
              caption = "Table 3. Presence/Ausence of all the isoform models")

# -------------------- 
# TABLE 3: SIRV metrics

if (lrgasp == TRUE){
  t3 <- DT::datatable(sirv.metrics,
                      extensions = "Buttons",
                      options = list(
                        pageLength = 15,
                        dom = 'Bfrtip',
                        buttons = c('copy', 'csv', 'pdf', 'print')
                      ),
                      escape = FALSE,
                      caption = "Table 4. SIRV metrics")
  t4 <- DT::datatable(fsm.metrics,
                      extensions = "Buttons",
                      options = list(
                        pageLength = 15,
                        dom = 'Bfrtip',
                        buttons = c('copy', 'csv', 'pdf', 'print')
                      ),
                      escape = FALSE,
                      caption = "Table 5. FSM metrics")
  t5 <- DT::datatable(ism.metrics,
                      extensions = "Buttons",
                      options = list(
                        pageLength = 15,
                        dom = 'Bfrtip',
                        buttons = c('copy', 'csv', 'pdf', 'print')
                      ),
                      escape = FALSE,
                      caption = "Table 6. ISM metrics")
  t6 <- DT::datatable(nic.metrics,
                      extensions = "Buttons",
                      options = list(
                        pageLength = 15,
                        dom = 'Bfrtip',
                        buttons = c('copy', 'csv', 'pdf', 'print')
                      ),
                      escape = FALSE,
                      caption = "Table 7. NIC metrics")
  t7 <- DT::datatable(nnc.metrics,
                      extensions = "Buttons",
                      options = list(
                        pageLength = 15,
                        dom = 'Bfrtip',
                        buttons = c('copy', 'csv', 'pdf', 'print')
                      ),
                      escape = FALSE,
                      caption = "Table 8. NNC metrics")
  
}


# -------------------- 
# PLOT 1: gene characterization
# PLOT 1.1: isoforms per gene

p1.1 <- ggplot(isoPerGene, aes(fill=number, y=countpergene, x=sample)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = myPalette) + mytheme 

# PLOT 1.2: exon structure

p1.2 <- ggplot(exonstructure, aes(fill=category, y=exonstructure, x=sample)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = myPalette) + mytheme

# PLOT 2: structural category distribution

p2 <- ggplot(df_SC, aes(fill=SC, y=value, x=ID)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = myPalette) + mytheme

# PLOT 3: splice junction distribution for each SC

#!!!!!!!!!!! ADD THESE PLOTS

# PLOT 4: distance to TSS

p4.1 <- ggplot(dist.list[[1]], aes(x=sample, y=log2(dist))) + 
  geom_boxplot(aes(col = sample)) +
  scale_fill_manual(values = myPalette) + mytheme

p4.2 <- ggplot(dist.list[[2]], aes(x=sample, y=log2(dist))) + 
  geom_boxplot(aes(col = sample)) +
  scale_fill_manual(values = myPalette) + mytheme

# PLOT5: distance to TTS

p5.1 <- ggplot(dist.list[[3]], aes(x=sample, y=log2(dist))) + 
  geom_boxplot(aes(col = sample)) +
  scale_fill_manual(values = myPalette) + mytheme

p5.2 <- ggplot(dist.list[[4]], aes(x=sample, y=log2(dist))) + 
  geom_boxplot(aes(col = sample)) +
  scale_fill_manual(values = myPalette) + mytheme

# PLOT 6: distance to CAGE peak

p6.1 <- ggplot(dist.list[[5]], aes(x=sample, y=log2(dist))) + 
  geom_boxplot(aes(col = sample)) +
  scale_fill_manual(values = myPalette) + mytheme

p6.2 <- ggplot(dist.list[[6]], aes(x=sample, y=log2(dist))) + 
  geom_boxplot(aes(col = sample)) +
  scale_fill_manual(values = myPalette) + mytheme

# PLOT 7: bad quality features
# PLOT 7.1: RT-switching

p7.1.1 <- ggplot(FSM.RT, aes(x=sample, y=FSM, fill=sample)) + geom_bar(stat="identity") +
  scale_fill_manual(values = myPalette) + mytheme
p7.1.2 <- ggplot(NIC.RT, aes(x=sample, y=NIC)) + geom_bar(color="blue", fill=rgb(0.1,0.4,0.5,0.7),stat="identity") +
  scale_fill_manual(values = myPalette) + mytheme
p7.1.3 <- ggplot(NNC.RT, aes(x=sample, y=NNC)) + geom_bar(color="blue", fill=rgb(0.1,0.4,0.5,0.7),stat="identity") +
  scale_fill_manual(values = myPalette) + mytheme

# p8: good quality features

# PLOT 9: Venn diagrams

# To not generate the .log files of VennDiagram
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

#p9 <- venn.diagram(l, filename = NULL, fill = myPalette[1:length(l)])

# PLOT 10: UpSet plot

p10 <-
  UpSetR::upset(
    UpSetR::fromList(l),
    order.by = "freq",
    mainbar.y.label = "Isoform Intersections",
    sets.x.label = "Isoforms Per Sample",
    nintersects = 20
    #sets.bar.color = myPalette[1:length(l)]
  )

# PLOT 11: Venn diagrams for SC

#p11 <- list()
#contador <- 1
#for (i in 1:length(res[[2]])) {
#  a <- res[[2]][res[[2]]$structural_category == str_cat[i], names(res[[2]])[3:length(names(res[[2]]))]]
#  l <- list()
#  for (j in 1:ncol(a)) {
#    l[[j]] <- na.omit(a[,j])
#  }
#  names(l) <- names(a)
#  p11[[contador]] <- venn.diagram(l, filename = NULL, fill = myPalette[1:length(l)])
#  contador <- contador + 1
#}

# PLOT 12: UpSet plots for SC

p12 <- list()
#for (i in 1:length(res$comparison)) {
#  a <- res$comparison[res$comparison$structural_category == str_cat[i], names(res$comparison)[3:length(names(res$comparison))]]
#  l <- list()
#  for (j in 1:ncol(a)) {
#    l[[j]] <- na.omit(a[,j])
#  }
#  names(l) <- names(a)
  
#  p12[[i]] <-
#    UpSetR::upset(
#      UpSetR::fromList(l),
#      order.by = "freq",
#      mainbar.y.label = "Isoform Intersections",
#      sets.x.label = "Isoforms Per Sample",
#      nintersects = 20
#     #sets.bar.color = myPalette[1:length(l)]
#    )
#}

if (TSS_TTS_coord == TRUE) {
  
  # PLOT 13: TSS standard deviation per pipeline
  
  a <- bind_rows(res$classifications, .id = "experiment_id")
  p13 <- ggplot(a, aes(log2(a[,"sdTSS"]))) +
    geom_density(aes(col = experiment_id)) + xlab(paste0("log2(",colnames(a)[8],")")) +
    scale_fill_manual(values = myPalette) + mytheme
  
  # PLOT 14: TTS standard deviation per pipeline
  
  p14 <- ggplot(a, aes(log2(a[,"sdTTS"]))) +
    geom_density(aes(col = experiment_id)) + xlab(paste0("log2(",colnames(a)[9],")")) +
    scale_fill_manual(values = myPalette) + mytheme
  

  # PLOT 15: TSS and TTS SD UJC
    p15.1 <-  ggplot(res$iso_metrics, aes(log2(res$iso_metrics$SD.TSS))) +
        geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
        mytheme + ggtitle(names(res$iso_metrics)[2]) 
    
    p15.2 <-  ggplot(res$iso_metrics, aes(log2(res$iso_metrics$SD.TTS))) +
      geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9) +
      mytheme + ggtitle(names(res$iso_metrics)[3]) 
    
    # Iso metrics
    # TSS
    
    p16.1 <- ggplot(res$iso_metrics, aes(x=res$iso_metrics$SD.TSS, y=res$iso_metrics$exons)) +
      geom_point(alpha=1/10) +
      geom_smooth(method="loess", color="red", se=FALSE) + xlab("SD TSS") + ylab("N exons") + ggtitle("UJC SD TSS vs Nº exons") + mytheme
    
    p16.2 <- ggplot(res$iso_metrics, aes(x=res$iso_metrics$SD.TSS, y=res$iso_metrics$length)) +
      geom_point(alpha=1/10) +
      geom_smooth(method="loess", color="red", se=FALSE)  + xlab("SD TSS") + ylab("N Length") + ggtitle("UJC SD TSS vs Length") + mytheme
    
    p16.3 <- ggplot(res$iso_metrics, aes(x=res$iso_metrics$SD.TSS, y=res$iso_metrics$FL_cpm)) +
      geom_point(alpha=1/10) +
      geom_smooth(method="loess", color="red", se=FALSE)  + xlab("SD TSS") + ylab("FL CPM") + ggtitle("UJC SD TSS vs CPM") + mytheme
    
    p16.4 <- ggplot(res$iso_metrics, aes(x=res$iso_metrics$length, y=res$iso_metrics$FL_cpm)) +
      geom_point(alpha=1/10) +
      geom_smooth(method="loess", color="red", se=FALSE)  + xlab("Length") + ylab("FL CPM") + ggtitle("UJC Length vs CPM") + mytheme
    
    # TTS
    p16.5 <- ggplot(res$iso_metrics, aes(x=res$iso_metrics$SD.TTS, y=res$iso_metrics$exons)) +
      geom_point(alpha=1/10) +
      geom_smooth(method="loess", color="red", se=FALSE) + xlab("SD TTS") + ylab("N exons") + ggtitle("UJC SD TTS vs Nº exons") + mytheme
    
    p16.6 <- ggplot(res$iso_metrics, aes(x=res$iso_metrics$SD.TTS, y=res$iso_metrics$length)) +
      geom_point(alpha=1/10) +
      geom_smooth(method="loess", color="red", se=FALSE)  + xlab("SD TTS") + ylab("N Length") + ggtitle("UJC SD TTS vs Length") + mytheme
    
    p16.7 <- ggplot(res$iso_metrics, aes(x=res$iso_metrics$SD.TTS, y=res$iso_metrics$FL_cpm)) +
      geom_point(alpha=1/10) +
      geom_smooth(method="loess", color="red", se=FALSE)  + xlab("SD TTS") + ylab("FL CPM") + ggtitle("UJC SD TTS vs CPM") + mytheme
    
    p16.8 <- ggplot(res$iso_metrics, aes(x=res$iso_metrics$length, y=res$iso_metrics$FL_cpm)) +
      geom_point(alpha=1/10) +
      geom_smooth(method="loess", color="red", se=FALSE)  + xlab("Length") + ylab("FL CPM") + ggtitle("UJC Length vs CPM") + mytheme
    
    # Isometris Log2
    
    p16.9 <- ggplot(res$iso_metrics, aes(x=log2(res$iso_metrics$SD.TSS), y=res$iso_metrics$exons)) +
      geom_point(alpha=1/10) +
      geom_smooth(method="loess", color="red", se=FALSE) + xlab("log 2 SD TSS") + ylab("N exons") + ggtitle("log2 SD TSS vs Nº exons") + mytheme
    
    p16.10 <- ggplot(res$iso_metrics, aes(x=log2(res$iso_metrics$SD.TSS), y=log2(res$iso_metrics$length))) +
      geom_point(alpha=1/10) +
      geom_smooth(method="loess", color="red", se=FALSE)  + xlab("log 2 SD TSS") + ylab("log 2 N Length") + ggtitle("log2 SD TSS vs Length") + mytheme
    
    p16.11 <- ggplot(res$iso_metrics, aes(x=log2(res$iso_metrics$SD.TSS), y=log2(res$iso_metrics$FL_cpm))) +
      geom_point(alpha=1/10) +
      geom_smooth(method="loess", color="red", se=FALSE)  + xlab("log 2 SD TSS") + ylab("log 2 FL CPM") + ggtitle("log2 SD TSS vs CPM") + mytheme
    
    p16.12 <- ggplot(res$iso_metrics, aes(x=log2(res$iso_metrics$length), y=log2(res$iso_metrics$FL_cpm))) +
      geom_point(alpha=1/10) +
      geom_smooth(method="loess", color="red", se=FALSE)  + xlab("log 2 Length") + ylab("log 2 FL CPM") + ggtitle("log2 Length vs CPM") + mytheme
    
    
    p16.13 <- ggplot(res$iso_metrics, aes(x=log2(res$iso_metrics$SD.TTS), y=res$iso_metrics$exons)) +
      geom_point(alpha=1/10) +
      geom_smooth(method="loess", color="red", se=FALSE) + xlab("log 2 SD TTS") + ylab("N exons") + ggtitle("log2 SD TTS vs Nº exons") + mytheme
    
    p16.14 <- ggplot(res$iso_metrics, aes(x=log2(res$iso_metrics$SD.TTS), y=log2(res$iso_metrics$length))) +
      geom_point(alpha=1/10) +
      geom_smooth(method="loess", color="red", se=FALSE)  + xlab("log 2 SD TTS") + ylab("log 2 N Length") + ggtitle("log2 SD TTS vs Length") + mytheme
    
    p16.15 <- ggplot(res$iso_metrics, aes(x=log2(res$iso_metrics$SD.TTS), y=log2(res$iso_metrics$FL_cpm))) +
      geom_point(alpha=1/10) +
      geom_smooth(method="loess", color="red", se=FALSE)  + xlab("log 2 SD TTS") + ylab("log 2 FL CPM") + ggtitle("log2 SD TTS vs CPM") + mytheme
    
    p16.16 <- ggplot(res$iso_metrics, aes(x=log2(res$iso_metrics$length), y=log2(res$iso_metrics$FL_cpm))) +
      geom_point(alpha=1/10) +
      geom_smooth(method="loess", color="red", se=FALSE)  + xlab("log 2 Length") + ylab("log 2 FL CPM") + ggtitle("log2 Length vs CPM") + mytheme
    
}

##*****
# Gene metrics

p17.1 <- ggplot(res$gene_metrics, aes(x=res$gene_metrics$length, y=res$gene_metrics$FL_cpm)) +
  geom_point(alpha=1/10) +
  geom_smooth(method="loess", color="red", se=FALSE) + xlab("Length") + ylab("FL CPM") + ggtitle("Gene Length vs CPM") + mytheme

p17.2 <- ggplot(res$gene_metrics, aes(x=res$gene_metrics$length, y=res$gene_metrics$N_UJC)) +
  geom_point(alpha=1/10) +
  geom_smooth(method="loess", color="red", se=FALSE) + xlab("Length") + ylab("N UJC") + ggtitle("Gene Length vs N_UJC") + mytheme

p17.3 <- ggplot(res$gene_metrics, aes(x=res$gene_metrics$length, y=res$gene_metrics$jaccard)) +
  geom_point(alpha=1/10) +
  geom_smooth(method="loess", color="red", se=FALSE) + xlab("Length") + ylab("Jaccard index")  + ggtitle("Gene Length vs Jaccard index") + mytheme

p17.4 <- ggplot(res$gene_metrics, aes(x=res$gene_metrics$N_UJC, y=res$gene_metrics$jaccard)) +
  geom_point(alpha=1/10) +
  geom_smooth(method="loess", color="red", se=FALSE) + xlab("N_UJC") + ylab("Jaccard index") + ggtitle("Gene N_UJC vs Jaccard index") + mytheme



p17.5 <- ggplot(res$gene_metrics, aes(x=log2(res$gene_metrics$length), y=log2(res$gene_metrics$FL_cpm))) +
  geom_point(alpha=1/10) +
  geom_smooth(method="loess", color="red", se=FALSE) + xlab("log 2 Length") + ylab("log 2 FL CPM") + ggtitle("log2 Gene Length vs CPM") + mytheme

p17.6 <- ggplot(res$gene_metrics, aes(x=log2(res$gene_metrics$length), y=res$gene_metrics$N_UJC)) +
  geom_point(alpha=1/10) +
  geom_smooth(method="loess", color="red", se=FALSE) + xlab("log 2 Length") + ylab("N UJC") + ggtitle("log 2 Gene Length vs N_UJC") + mytheme

p17.7 <- ggplot(res$gene_metrics, aes(x=log2(res$gene_metrics$length), y=res$gene_metrics$jaccard)) +
  geom_point(alpha=1/10) +
  geom_smooth(method="loess", color="red", se=FALSE) + xlab("log 2 Length") + ylab("Jaccard index")  + ggtitle("log 2 Gene Length vs Jaccard index") + mytheme


# Gene metrics log2

# -------------------- Output report

print("Generating the HTML report...")
#Sys.setenv(RSTUDIO_PANDOC = "/usr/lib/rstudio/bin/pandoc")
rmarkdown::render(
  input = paste(getwd(), "SQANTI3_comparison_report.Rmd", sep = "/"),
  output_dir = output_directory,
  output_file = paste0(output_name, ".html")
)

# -------------------- Output csv

print("Writting CSV file with P/A table of UJC...")
write.csv(res$comparisonPA, paste0(output_directory, "/", output_name, ".csv"))

print("DONE\nExecution finished")

