#####################################
##### SQANTI3 report generation ######
#####################################



### Author: Lorena de la Fuente, Elizabeth Tseng & Francisco J Pardo-Palacios
### Last Modified: 09/23/2020 by etseng@pacb.com

#********************** Taking arguments from python script

args <- commandArgs(trailingOnly = TRUE)
class.file <- args[1]
junc.file <- args[2]
utilities.path <- args[4]
saturation.curves <- args[5]
report.format <- args[6]

if (length(args) < 6) {
  stop("Incorrect number of arguments! Script usage is: [classification file] [junction file] [utilities directory path] [True/False for saturation curves] [pdf|html|both]. Abort!")
}

if (!(saturation.curves %in% c('True', 'False'))) {
  stop("Saturation curve argument needs to be 'True' or 'False'. Abort!")
}

if (!(report.format %in% c('pdf', 'html', 'both'))) {
  stop("Report format needs to be: pdf, html, or both. Abort!")
}

source(paste(utilities.path, "/report_qc/generatePDFreport.R", sep = "/"))

if (saturation.curves=='True'){
  source(paste(utilities.path, "/report_qc/saturation/plot_saturation.R", sep = "/"))
  source(paste(utilities.path, "/report_qc/saturation/LR_saturation.R", sep = "/"))
  source(paste(utilities.path, "/report_qc/saturation/data_prep_saturation.R", sep = "/"))
}

report.prefix <- strsplit(class.file, "_classification.txt")[[1]][1];
output_directory <- dirname(class.file)
output_name <- basename(report.prefix)
pdf.report.file <- paste0(report.prefix, "_SQANTI3_report.pdf");
html.report.file <- paste0(output_name, "_SQANTI3_report.html");
class.file2 <- paste(report.prefix, "_classification_TPM.txt", sep='');


#********************** Packages 

library(ggplot2)
library(ggplotify)
library(scales)
library(reshape)
library(grid)
library(gridExtra)
library(NOISeq)
library(rmarkdown)
library(htmltools)
library(DT)
library(plyr)
library(plotly)
library(dplyr)

# ***********************
# ***********************
# PLOTS INDEX


# PLOT length of isoforms
# p.length.all: length of all isoforms, regardless of category
# p.length.cat: length of isoforms, by category
# p.length.exon: length of isoforms, mono- vs mult-exon
# (optional) p.length.all.sample
# (optional) p.length.exon.sample


# p0: Splicing complexity (X) Isoforms per Gene (Y) Number of Genes
# p1: Distribution of Isoform Classification
# p2: Distribution of Ref Lengths (FSM ISM only)
# p3: Distribution of Ref Exons   (FSM ISM only)
# p4: Distribution of Isoform Lengths, by Classification
# p5: Distribution of Exon Counts, by Classification
# p6: Distribution of Mono- vs Multi-Exons, Novel vs Annotated Genes
# p7: Splicing complexity, Novel vs Annotated Genes
# p8: Transcript Expression (SR) by Structural Category { log2(TPM+1) }
# p9: Transcript Expression (FL) by Structural Category { log2(TPM+1) }
# p10: Gene Expression (SR), Annotated vs Novel genes { log2(Gene_TPM+1) }
# p11: Gene Expression (SR), Annotated vs Novel genes { log2(Gene_TPM+1) }
# p13: Gene Expression level in NNC/FSM containing genes
# p13.c: Gene Expression level in NNC/FSM containing genes

# p.classByLen.a: Structural categories with increasing transcript length, absolute
# p.classByLen.b: Structural categories with increasing transcript length, normalized
#
# p21.a: Distance to polyA site for FSM, absolute
# p21.b: Distance to polyA site for FSM, percentage
# p21.dist3.ISM.a: Distance to polyA site for ISM, absolute
# p21.dist3.ISM.b: Distance to polyA site for ISM, percentage
# p22.a: Distance to start site for FSM, absolute
# p22.b: Distance to start site for FSM, percentage

# p23.a: Splice Junctions by Classification (known canonical, known non-canonical, novel canonical, novel non-canonical)
# p23.b: Splice Junctions by Classification (canonical vs non-canonical)
#
# p28.a: Good Quality Control Attributes Across Structural Categories (annot, SJ, coverage)
# p28.aa: Good Quality Control Attributes Across Structural Categories (polyA, Cage)
# p28.a.SJ: Percentage of  All Canonical Junctions
# p28.a.Cov: Percentage of Splice Junctions With Short Reads Coverage
# p28.a.Cage: Percentage of Cage Support
# p28.a.polyA : Percentage of PolyA Support
# p28.a.annot : Percentage of Annotation Support
# p29.a: Splice Junction, % of RT switching, all junctions
# p29.b: Splice Junction, % of RT switching, unique junctions
#
# p30: intra-priming, by Classification
# p31: intra-priming, Mono- vs Multi-Exons
# p32: intra-priming, Coding vs Non-Coding


# ***********************

#********************** Reading Information

########## Classification information

data.class = read.table(class.file, header=T, as.is=T, sep="\t")
rownames(data.class) <- data.class$isoform

xaxislevelsF1 <- c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog", "genic","antisense","fusion","intergenic","genic_intron");
xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")
subc.levels=c("alternative_3end",'alternative_3end5end', "alternative_5end","reference_match", "3prime_fragment","internal_fragment", "5prime_fragment","combination_of_known_junctions", "combination_of_known_splicesites", "intron_retention","no_combination_of_known_junctions", "mono-exon_by_intron_retention", "at_least_one_novel_splicesite", "mono-exon", "multi-exon")
subc.labels=c("Alternative 3'end", "Alternative 3'5'end", "Alternative 5'end", "Reference match", "3' fragment", "Internal fragment", "5' fragment", "Comb. of annot. junctions", "Comb. of annot. splice sites", "Intron retention", "Not comb. of annot. junctions", "Mono-exon by intron ret.", "At least 1 annot. don./accept.", "Mono-exon", "Multi-exon")
coding.levels=c("coding", "non_coding")
coding.labels=c("Coding", "Non coding")


data.class$structural_category = factor(data.class$structural_category,
                                        labels = xaxislabelsF1,
                                        levels = xaxislevelsF1,
                                        ordered=TRUE)
data.class$subcategory = factor(data.class$subcategory,
                                        labels = subc.labels,
                                        levels = subc.levels,
                                        ordered=TRUE)
data.class$coding = factor(data.class$coding,
                                labels = coding.labels,
                                levels = coding.levels,
                                ordered=TRUE)
legendLabelF1 <- levels(as.factor(data.class$coding));

####### SRTM and SNTM functions

STM_function <- function(x){
  five=FALSE
  three=FALSE
  sj=FALSE
  ev_sj <- !is.na(x["min_cov"]) & as.numeric(x["exons"])>1
  ref_TSS <- FALSE
  ref_TTS <- FALSE
  
  if (!is.na(x["diff_to_gene_TSS"])){
    if (abs(as.numeric(x["diff_to_gene_TSS"]))<=50){
      ref_TSS <- TRUE
    }
  }
  
  if (!is.na(x["diff_to_gene_TTS"])){
    if (abs(as.numeric(x["diff_to_gene_TTS"]))<=50){
      ref_TTS <- TRUE
    }
  }
  
  w_cage <- !is.na(x["within_CAGE_peak"]) & x["within_CAGE_peak"]=="True"
  if ( ref_TSS | w_cage  ){
    five=TRUE
  }
  w_polya <- !is.na(x["within_polyA_site"]) & x["within_polyA_site"]=="True"
  if (ref_TTS | w_polya | !is.na(x["polyA_motif"])){
    three=TRUE
  }
  if (x["structural_category"]=="FSM" | x["structural_category"]=="ISM"){
    sj=TRUE
  }else if ( ev_sj ){
    if (as.numeric(x["min_cov"])>0){
      sj=TRUE
    }
  }
  
  if (five & three & sj){
    return("Fully supported")
  }else{
    return("Not fully supported")
  }
}

data.class$STM <- apply(data.class,1, STM_function)
################################

data.FSMISM <- subset(data.class, structural_category %in% c('FSM', 'ISM'))
data.NICNNC <- subset(data.class, structural_category %in% c("NIC", "NNC"))
data.other <- subset(data.class, structural_category %in% c("Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron"))
data.FSM <- subset(data.class, (structural_category=="FSM" & exons>1))
data.ISM <- subset(data.class, (structural_category=="ISM" & exons>1))
data.NNC <- subset(data.class, (structural_category=="NNC" & exons>1))
data.NIC <- subset(data.class, (structural_category=="NIC" & exons>1))
data.GenicGenomic <- subset(data.class, (structural_category=="Genic\nGenomic" & exons>1 ))
data.Antisense <- subset(data.class, (structural_category=="Antisense" & exons>1))
data.Fusion <- subset(data.class, (structural_category=="Fusion" & exons>1))
data.Intergenic <- subset(data.class, (structural_category=="Intergenic" & exons>1))
data.GenicIntron <- subset(data.class, (structural_category=="Genic\nIntron" & exons>1))

# subcategories data sets
#"FSM"
data.alt3end <- subset(data.FSM, (subcategory=="Alternative 3'end"))
data.alt35end <- subset(data.FSM, (subcategory=="Alternative 3'5'end"))
data.alt5end <- subset(data.FSM, (subcategory=="Alternative 5'end"))
data.refmatch <- subset(data.FSM, (subcategory=="Reference match"))
#"ISM"
data.3fragment <- subset(data.ISM, (subcategory=="3' fragment"))
data.int_fragment <- subset(data.ISM, (subcategory=="Internal fragment"))
data.5fragment <- subset(data.ISM, (subcategory=="5' fragment"))
data.intron_ret_ISM <- subset(data.ISM, (subcategory=="Intron retention"))
#"NIC"
data.comb_annot_js_NIC <- subset(data.NIC, (subcategory=="Comb. of annot. junctions"))
data.comb_annot_ss_NIC <- subset(data.NIC, (subcategory=="Comb. of annot. splice sites"))
data.intron_ret_NIC <- subset(data.NIC, (subcategory=="Intron retention"))
data.mono_ex_intron_ret_NIC <- subset(data.NIC, (subcategory=="Mono-exon by intron ret."))
#"NNC"
data.comb_annot_js_NNC <- subset(data.NNC, (subcategory=="Comb. of annot. junctions"))
data.comb_annot_ss_NNC <- subset(data.NNC, (subcategory=="Comb. of annot. splice sites"))
data.intron_ret_NNC <- subset(data.NNC, (subcategory=="Intron retention"))
data.mono_ex_intron_ret_NNC <- subset(data.NNC, (subcategory=="Mono-exon by intron ret."))
data.one_don_acc <- subset(data.NNC, (subcategory=="At least 1 annot. don./accept."))


subcategories.FSM <- list(data.alt3end, data.alt35end, data.alt5end, data.refmatch)
subcategories.ISM <- list(data.3fragment, data.int_fragment, data.5fragment, data.intron_ret_ISM)
subcategories.NIC <- list(data.comb_annot_js_NIC, data.comb_annot_ss_NIC, data.intron_ret_NIC, data.mono_ex_intron_ret_NIC)
subcategories.NNC <- list(data.comb_annot_js_NNC, data.comb_annot_ss_NNC, data.intron_ret_NNC, data.mono_ex_intron_ret_NNC, data.one_don_acc)

########### Junction information
data.junction <- read.table(junc.file, header=T, as.is=T, sep="\t")

# make a unique identifier using chrom_strand_start_end
data.junction$junctionLabel = with(data.junction, paste(chrom, strand, genomic_start_coord, genomic_end_coord, sep="_"))

data.junction$SJ_type <- with(data.junction, paste(junction_category,canonical,"SJ", sep="_"))
data.junction$SJ_type <- factor(data.junction$SJ_type, levels=c("known_canonical_SJ", "known_non_canonical_SJ", "novel_canonical_SJ", "novel_non_canonical_SJ"),
                                labels=c("Known\ncanonical ", "Known\nNon-canonical ", "Novel\ncanonical ", "Novel\nNon-canonical "), order=T)

data.junction$structural_category = data.class[data.junction$isoform, "structural_category"]

uniqJunc <- unique(data.junction[,c("junctionLabel", "SJ_type", "total_coverage_unique")]);
uniqJuncRTS <- unique(data.junction[,c("junctionLabel","SJ_type", "RTS_junction")]);

########## Generating plots

#*** Global plot parameters

myPalette = c("#6BAED6","#FC8D59","#78C679","#EE6A50","#969696","#66C2A4", "goldenrod1", "darksalmon", "#41B6C4","tomato3", "#FE9929")
subcat.palette = c("Alternative 3'end"='#02314d',
                   "Alternative 3'5'end"='#0e5a87',
                   "Alternative 5'end"='#7ccdfc',
                   'Reference match'='#c4e1f2',
                   "3' fragment"='#c4531d',
                   "Internal fragment"='#e37744',  
                   "5' fragment"='#e0936e', 
                   "Comb. of annot. junctions"='#014d02',
                   "Comb. of annot. splice sites"='#379637',  
                   "Intron retention"='#81eb82', 
                   "Not comb. of annot. junctions"='#6ec091',
                   "Mono-exon by intron ret."='#4aaa72',
                   "At least 1 annot. don./accept."='#32734d',
                   "Mono-exon"='#cec2d2',
                   "Multi-exon"='#876a91')



cat.palette = c("FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                "NNC"="#EE6A50", "Genic\nGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                "Intergenic" = "darksalmon", "Genic\nIntron"="#41B6C4")
 

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


# Create a new attribute called "novelGene"

data.class$novelGene <- "Annotated Genes"
data.class[grep("novelGene", data.class$associated_gene), "novelGene"] <- "Novel Genes"
data.class$novelGene = factor(data.class$novelGene,
                              levels = c("Novel Genes","Annotated Genes"),
                              ordered=TRUE)

# Create a new attribute called "exonCat"

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

# ----------------------------------------------------------
# Make "isoPerGene" which is aggregated information by gene
#  $associatedGene - either the ref gene name or novelGene_<index>
#  $novelGene      - either "Novel Genes" or "Annotated Genes"
#  $FSM_class      - "A", "B", or "C"
#  $geneExp        - gene expression info
#  $nIso           - number of isoforms associated with this gene
#  $nIsoCat        - splicing complexity based on number of isoforms
# ----------------------------------------------------------
if (!all(is.na(data.class$gene_exp))){
  isoPerGene = aggregate(data.class$isoform,
                         by = list("associatedGene" = data.class$associated_gene,
                                   "novelGene" = data.class$novelGene,
                                   "FSM_class" = data.class$FSM_class,
                                   "geneExp"=data.class$gene_exp),
                         length)
} else {
  isoPerGene = aggregate(data.class$isoform,
                         by = list("associatedGene" = data.class$associated_gene,
                                   "novelGene" = data.class$novelGene,
                                   "FSM_class" = data.class$FSM_class),
                         length)
}
# assign the last column with the colname "nIso" (number of isoforms)
colnames(isoPerGene)[ncol(isoPerGene)] <- "nIso"


isoPerGene$FSM_class2 = factor(isoPerGene$FSM_class,
                               levels = c("A", "B", "C"),
                               labels = c("MonoIsoform Gene", "MultiIsoform Genes\nwithout expression\nof a FSM", "MultiIsoform Genes\nexpressing at least\none FSM"),
                               ordered=TRUE)

isoPerGene$novelGene = factor(isoPerGene$novelGene,
                              levels = c("Annotated Genes", "Novel Genes"),
                              ordered=TRUE)

max_iso_per_gene <- max(isoPerGene$nIso)
if (max_iso_per_gene >= 6) {
    isoPerGene$nIsoCat <- cut(isoPerGene$nIso, breaks = c(0,1,3,5,max_iso_per_gene+1), labels = c("1", "2-3", "4-5", ">=6"));
} else if (max_iso_per_gene >= 5) {
    isoPerGene$nIsoCat <- cut(isoPerGene$nIso, breaks = c(0,1,3,5), labels = c("1", "2-3", "4-5"));
} else if (max_iso_per_gene >= 3) {
    isoPerGene$nIsoCat <- cut(isoPerGene$nIso, breaks = c(0,1,3), labels = c("1", "2-3"));
} else {
    isoPerGene$nIsoCat <- cut(isoPerGene$nIso, breaks = c(0,1), labels = c("1"));
}

# see if there are multple FL samples
FL_multisample_indices <- which(substring(colnames(data.class), 1,3)=="FL.")

if (any(grep('^FL\\.', names(data.class)))){
  data.class$FL= rowSums(data.class[,grep('^FL\\.', names(data.class))])
}

if (!all(is.na(data.class$FL))){
  total_fl <- sum(data.class$FL, na.rm=T)
  #data.class$FL_TPM <- round(data.class$FL*(10**6)/total_fl)
  data.class$FL_TPM <- data.class$FL*(10**6)/total_fl
}

if (!all(is.na(data.FSMISM$FL))){
  total_fl <- sum(data.FSMISM$FL, na.rm=T)
  #data.class$FL_TPM <- round(data.class$FL*(10**6)/total_fl)
  data.FSMISM$FL_TPM <- data.FSMISM$FL*(10**6)/total_fl
}

if (!all(is.na(data.NICNNC$FL))){
  total_fl <- sum(data.NICNNC$FL, na.rm=T)
  #data.class$FL_TPM <- round(data.class$FL*(10**6)/total_fl)
  data.NICNNC$FL_TPM <- data.NICNNC$FL*(10**6)/total_fl
}

if (!all(is.na(data.other$FL))){
  total_fl <- sum(data.other$FL, na.rm=T)
  #data.class$FL_TPM <- round(data.class$FL*(10**6)/total_fl)
  data.other$FL_TPM <- data.other$FL*(10**6)/total_fl
}

if (length(FL_multisample_indices)>0){
  FL_multisample_names <- substring(colnames(data.class)[FL_multisample_indices],4)
  FL_TPM_multisample_names <- c();

  for (i in 1:length(FL_multisample_indices)) {
    j <- FL_multisample_indices[i]
    name <- paste("FL_TPM", FL_multisample_names[i], sep='.')
    name2 <- paste(name, "_log10", sep='')
    FL_TPM_multisample_names <- c(FL_TPM_multisample_names, name)
    total_fl <- sum(data.class[j])
    data.class[,name] <- data.class[j]*(10**6)/total_fl
    data.class[,name2] <- log10(data.class[j]*(10**6)/total_fl + 1)
  }
  data.class$structural_category = factor(data.class$structural_category,
                                        labels = xaxislevelsF1,
                                        levels = xaxislabelsF1,
                                        ordered=TRUE)
  write.table(data.class, class.file2, quote=F, sep = "\t", row.names = FALSE);
  data.class$structural_category = factor(data.class$structural_category,
                                        labels = xaxislabelsF1,
                                        levels = xaxislevelsF1,
                                        ordered=TRUE)
}


########################################
######### LENGTH PLOTS  ################
########################################

p.length.all <- ggplot(data.class, aes(x=length)) +
  geom_histogram(binwidth=100) +
  labs(x="Transcript length", y="Count", title="All Transcript Lengths Distribution") +
  theme(legend.position="top") +
  scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
  mytheme

if (length(FL_multisample_indices)>0) {  # has multiple samples
  df.length_by_sample <- data.frame();
  for (i in 1:length(FL_multisample_indices)) {
    j <- FL_multisample_indices[i];
    df_new <- data.class[data.class[j]>0,c("isoform", "length", "exonCat")];
    df_new$sample <- FL_multisample_names[i];
    df.length_by_sample <- rbind(df.length_by_sample, df_new);
  }

  p.length.all.sample <- ggplot(df.length_by_sample, aes(x=length, color=sample)) +
    geom_freqpoly(binwidth=100) +
    scale_color_manual(values = cat.palette, name="Structural Category") +
    labs(x="Transcript length", y="Count", title="Transcript Lengths by Sample") +
    theme(legend.position="bottom") +
    scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
    mytheme

  p.length.exon.sample <- ggplot(df.length_by_sample, aes(x=length, color=sample, lty=exonCat)) +
    geom_freqpoly(binwidth=100) +
    labs(x="Transcript length", y="Count", title="Mono- vs Multi-Exons Transcript Lengths by Sample") +
    theme(legend.position="bottom") +
    scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
    mytheme
}

p.length.cat <- ggplot(data.class, aes(x=length, color=structural_category)) +
  geom_freqpoly(binwidth=100, size=1) +
  scale_color_manual(values = cat.palette, name="Structural Category") +
  labs(x="Transcript length", y="Count", title="Transcript Lengths Distribution by Structural Category") +
  scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
  mytheme+
  theme(legend.position="bottom", legend.title=element_blank())

p.length.exon <- ggplot(data.class, aes(x=length, color=exonCat)) +
  geom_freqpoly(binwidth=100, size=1) +
  labs(x="Transcript length", y="Count", title="Mono- vs Multi- Exon Transcript Lengths Distribution") +
  scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
  mytheme +
  theme(legend.position="bottom", legend.title=element_blank()) 


#**** PLOT 0: Distribution of Number of Isoforms per Gene

p0 <- ggplot(isoPerGene, aes(x=nIsoCat, fill=nIsoCat)) +
  geom_bar(stat="count", aes(y= (..count..)/sum(..count..)*100), color="black", size=0.3, width=0.7) +
  guides(fill="none") +
  scale_y_continuous(labels = function(x) format(x), expand = c(0,0)) +
  scale_fill_manual(values = myPalette[c(2:5)]) +
  labs(x ="Isoforms per gene", title="Number of Isoforms per Gene\n\n\n", y = "Genes, %") +
  mytheme

#**** PLOT 1: Structural Classification

p1 <- ggplot(data=data.class, aes(x=structural_category)) +
  geom_bar(aes(y = (..count..)/sum(..count..)*100, alpha=coding, fill=structural_category), color="black", size=0.3, width=0.7) +
  #geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = "count", vjust = -0.25)  +
  scale_x_discrete(drop=FALSE) +
  scale_alpha_manual(values=c(1,0.3),
                     name = "Coding prediction",
                     labels = legendLabelF1)+
  xlab("") +
  ylab("Transcripts, %") +
  mytheme +
  geom_blank(aes(y=((..count..)/sum(..count..))), stat = "count") +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = cat.palette, guide='none') +
  ggtitle("Isoform Distribution Across Structural Categories\n\n" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12)) +
  scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))

p1.s.titles = list("Isoform Distribution Across FSM\n\n",
                   "Isoform Distribution Across ISM\n\n",
                   "Isoform Distribution Across NNC\n\n",
                   "Isoform Distribution Across NIC\n\n",
                   "Isoform Distribution Across Genic Genomic\n\n",
                   "Isoform Distribution Across Antisense\n\n",
                   "Isoform Distribution Across Fusion\n\n",
                   "Isoform Distribution Across Intergenic\n\n",
                   "Isoform Distribution Across Genic Intron\n\n")

categories.list=list(data.FSM, data.ISM, data.NNC, data.NIC, data.GenicGenomic, data.Antisense, 
                     data.Fusion, data.Intergenic, data.GenicIntron)

p1.s.list = list()
for(i in 1:length(categories.list)){
  c<-data.frame(categories.list[i])
  if (!(dim(c))[1]==0){
    p1.s <- ggplot(data=c, aes(x=subcategory)) +
      geom_bar(aes(y = (..count..)/sum(..count..)*100, alpha=coding, fill=subcategory), color="black", size=0.3, width=0.7) +
      scale_x_discrete(drop=TRUE) +
      scale_alpha_manual(values=c(1,0.3), name = "Coding prediction", labels = legendLabelF1)+
      ylab("Transcripts, %") +
      mytheme +
      geom_blank(aes(y=((..count..)/sum(..count..))), stat = "count") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_manual(values = subcat.palette, guide='none') +
      ggtitle(p1.s.titles[i])+
      scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
      theme(axis.title.x=element_blank())  
    p1.s.list[[i]] = p1.s
  }
}

#**** PLOTS 2-3: refLength and refExons for ISM and FSM transcripts. Plot if any ISM or FSM transcript

if (nrow(data.FSMISM) > 0) {

  p2 <- ggplot(data=data.FSMISM, aes(x=structural_category, y=ref_length/1000, fill=structural_category)) +
    geom_boxplot(color="black", size=0.3, outlier.size = 0.2) + mytheme +
    scale_fill_manual(values = myPalette) +
    scale_x_discrete(drop = TRUE) +
    guides(fill="none") +
    xlab("") +
    ylab("Matched reference length, kb") +
    labs(title="Length Distribution of Matched Reference Transcripts\n\n\n",
         subtitle="Applicable Only to FSM and ISM Categories\n\n")

  p3 <- ggplot(data=data.FSMISM, aes(x=structural_category, y=ref_exons, fill=structural_category)) +
    geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
    scale_x_discrete(drop = TRUE) +
    xlab("") +
    ylab("Matched reference exon count") +
    scale_fill_manual(values = myPalette) +
    guides(fill="none") +
    mytheme +
    labs(title="Exon Count Distribution of Matched Reference Transcripts\n\n\n",
         subtitle="Applicable Only to FSM and ISM Categories\n\n")

}

#****  PLOT 4: Transcript lengths by category

p4 <- ggplot(data=data.class, aes(x=structural_category, y=length, fill=structural_category)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  scale_x_discrete(drop=FALSE) +
  ylab("Transcript Length (bp)") +
  scale_fill_manual(values = cat.palette) +
  guides(fill="none") +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) +
  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
  ggtitle("Transcript Lengths by Structural Classification\n\n" ) +
  theme(axis.title.x=element_blank())

p4.s1 <- ggplot(data=data.FSMISM, aes(x=structural_category, y=length, fill=subcategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  scale_x_discrete(drop=TRUE) +
  ylab("Transcript Length (bp)") +
  scale_fill_manual(values = subcat.palette, drop=TRUE) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.position="right", legend.title=element_blank()) +
  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
  ggtitle("Transcript Lengths by Subcategory\n\n" ) +
  theme(axis.title.x=element_blank())

p4.s2 <- ggplot(data=data.NICNNC, aes(x=structural_category, y=length, fill=subcategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  scale_x_discrete(drop=TRUE) +
  ylab("Transcript Length (bp)") +
  scale_fill_manual(values = subcat.palette, drop=TRUE) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.position="right", legend.title=element_blank())+
  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
  ggtitle("Transcript Lengths by Subcategory\n\n" ) +
  theme(axis.title.x=element_blank())

p4.s3 <- ggplot(data=data.other, aes(x=structural_category, y=length, fill=subcategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  scale_x_discrete(drop=TRUE) +
  ylab("Transcript Length (bp)") +
  scale_fill_manual(values = subcat.palette, drop=TRUE) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) +
  theme(legend.position="right", legend.title=element_blank())+
  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
  ggtitle("Transcript Lengths by Subcategory\n\n" ) +
  theme(axis.title.x=element_blank())


##**** PLOT 5: Exon counts by category

p5 <- ggplot(data=data.class, aes(x=structural_category, y=exons, fill=structural_category)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  ylab("Number of exons") +
  scale_x_discrete(drop=FALSE) +
  scale_fill_manual(values = cat.palette) +
  guides(fill="none") +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) +
  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
  ggtitle("Exon Counts by Structural Classification\n\n" ) +
  theme(axis.title.x=element_blank())

###Exon counts by subcategory
p5.s1 <- ggplot(data=data.FSMISM, aes(x=structural_category, y=exons, fill=subcategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) +
  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
  theme(legend.position="right", legend.title=element_blank()) +
  theme(axis.title.x=element_blank())+
  ylab("Number of exons") +
  scale_x_discrete(drop=TRUE) +
  scale_fill_manual(values = subcat.palette) +
  ggtitle("Exon Counts by Subcategory\n\n" )

p5.s2 <- ggplot(data=data.NICNNC, aes(x=structural_category, y=exons, fill=subcategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  ylab("Number of exons") +
  scale_x_discrete(drop=TRUE) +
  scale_fill_manual(values = subcat.palette) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) +
  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
  theme(legend.position="right", legend.title=element_blank())+
  ggtitle("Exon Counts by Subcategory\n\n" ) +
  theme(axis.title.x=element_blank())

p5.s3 <- ggplot(data=data.other, aes(x=structural_category, y=exons, fill=subcategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  ylab("Number of exons") +
  scale_x_discrete(drop=TRUE) +
  scale_fill_manual(values = subcat.palette) +
  mytheme  + theme(axis.text.x = element_text(angle = 45)) +
  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
  theme(legend.position="right", legend.title=element_blank())+
  ggtitle("Exon Counts by Subcategory\n\n" ) +
  theme(axis.title.x=element_blank())

##### STM plots

data.FSMISMNICNNC <- rbind(data.FSMISM, data.NICNNC )
pSTM <- ggplot(data=data.FSMISMNICNNC, aes(x=structural_category)) +
  geom_bar(aes(y = (..count..), alpha=STM, fill=structural_category), color="black", size=0.3, width=0.7) +
  scale_x_discrete(drop=TRUE) +
  scale_alpha_manual(values=c(1,0.3),
                     name = "Supported Transcript Model",
                     guide = "legend")+
  xlab("") +
  ylab("Transcripts, count") +
  mytheme +
  geom_blank(aes(y=(..count..)), stat = "count") +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = cat.palette, guide='none') +
  ggtitle("Isoform Distribution Across Structural Categories\n\n" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12)) +
  scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
  theme(legend.position = "right")

pSTM_perc <- ggplot(data=data.FSMISMNICNNC, aes(x=structural_category)) +
  geom_bar(aes(y = (..count..), alpha=STM, fill=structural_category), position="fill", color="black", size=0.3, width=0.7) +
  scale_x_discrete(drop=TRUE) +
  scale_alpha_manual(values=c(1,0.3),
                     name = "Supported Transcript Model",
                     guide = "legend")+
  xlab("") +
  ylab("Transcripts, %") +
  mytheme +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = cat.palette, guide='none') +
  ggtitle("Isoform Distribution Across Structural Categories\n\n" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12)) +
  scale_y_continuous(expand=expansion(mult = c(0,0.1)), labels = scales::percent) +
  theme(legend.position = "right")


pSTM.s1 <- ggplot(data=data.FSMISM, aes(x=subcategory)) +
  geom_bar(aes(y = (..count..), alpha=STM, fill=subcategory), color="black", size=0.3, width=0.7) +
  scale_x_discrete(drop=TRUE) +
  scale_alpha_manual(values=c(1,0.3),
                     name = "Supported Transcript Model",
                     guide = "legend")+
  xlab("") +
  ylab("Transcripts, count") +
  mytheme +
  facet_grid(.~ structural_category, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = subcat.palette, guide = "none") +
  ggtitle("Isoform Distribution Across Structural Subcategories\n\n",
          subtitle = "FSM and ISM" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x=element_text(size=10)) +
  scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
  theme(legend.position = "right")

pSTM_perc.s1 <- ggplot(data=data.FSMISM, aes(x=subcategory)) +
  geom_bar(aes(y = (..count..), alpha=STM, fill=subcategory), position="fill", color="black", size=0.3, width=0.7) +
  scale_x_discrete(drop=TRUE) +
  scale_alpha_manual(values=c(1,0.3),
                     name = "Supported Transcript Model",
                     guide = "legend")+
  xlab("") +
  ylab("Transcripts, %") +
  mytheme +
  facet_grid(.~ structural_category, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = subcat.palette, guide = "none") +
  ggtitle("Isoform Distribution Across Structural Subcategories\n\n",
          subtitle = "FSM and ISM" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(size=10)) +
  scale_y_continuous(expand=expansion(mult = c(0,0)), labels = scales::percent) +
  theme(legend.position = "right")

pSTM.s2 <- ggplot(data=data.NICNNC, aes(x=subcategory)) +
  geom_bar(aes(y = (..count..), alpha=STM, fill=subcategory), color="black", size=0.3, width=0.7) +
  scale_x_discrete(drop=TRUE) +
  scale_alpha_manual(values=c(1,0.3),
                     name = "Supported Transcript Model",
                     guide = "legend")+
  xlab("") +
  ylab("Transcripts, count") +
  mytheme +
  facet_grid(.~ structural_category, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = subcat.palette, guide='none') +
  ggtitle("Isoform Distribution Across Structural Subcategories\n\n",
          subtitle = "NIC and NNC") +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(size=10)) +
  scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
  theme(legend.position = "right")

pSTM_perc.s2 <- ggplot(data=data.NICNNC, aes(x=subcategory)) +
  geom_bar(aes(y = (..count..), alpha=STM, fill=subcategory), position="fill", color="black", size=0.3, width=0.7) +
  scale_x_discrete(drop=TRUE) +
  scale_alpha_manual(values=c(1,0.3),
                     name = "Supported Transcript Model",
                     guide = "legend")+
  xlab("") +
  ylab("Transcripts, %") +
  mytheme +
  facet_grid(.~ structural_category, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = subcat.palette, guide='none') +
  ggtitle("Isoform Distribution Across Structural Subcategories\n\n",
          subtitle = "NIC and NNC" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(size=10)) +
  scale_y_continuous(expand=expansion(mult = c(0,0)), labels = scales::percent) +
  theme(legend.position = "right")

##**** PLOT 6: Mono vs Multi-exon distribution for Known vs Novel Genes

p6 <- ggplot(data=data.class, aes(x=novelGene)) +
  geom_bar(position="fill",aes(y = (..count..)/sum(..count..), fill=exonCat), color="black", size=0.3, width=0.5) +
  scale_x_discrete(drop=FALSE) +
  scale_y_continuous(breaks=c(0.0,0.25,0.5,0.75,1),
                     labels=c("0","25","50","75","100"), expand = c(0,0)) +
  scale_fill_manual(name = "Transcript type",
                    values = myPalette[c(2:5)]) +
  ylab("Transcripts, %") +
  mytheme +
  theme(axis.title.x=element_blank()) +
  theme(legend.position="bottom") +
  ggtitle("Distribution of Mono- vs Multi-Exon Transcripts\n\n" )

# p7: Distribution of Number of Isoforms, separated by Novel vs Annotated Genes

p7 <- ggplot(data=isoPerGene, aes(x=novelGene)) +
  geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=nIsoCat), color="black", size=0.3, width=0.5) +
  scale_y_continuous(breaks=c(0.0,0.25,0.5,0.75,1),
                   labels=c("0","25","50","75","100"), expand = c(0,0)) +
  scale_x_discrete(drop=FALSE) +
  scale_fill_manual(name = "Isoforms Per Gene",
                    values = myPalette[c(2:5)]) +
  ylab("Genes, %") +
  xlab("Gene Type") +
  mytheme +
  theme(axis.title.x=element_blank()) +
  theme(legend.position="bottom") +
  guides(fill = guide_legend(keywidth = 0.9, keyheight = 0.9)) +
  labs(title="Number of Isoforms per Gene\n\n\n",
       subtitle="Known vs Novel Genes\n\n")


##**** PLOT  absolute and normalized % of different categories with increasing transcript length
# requires use of dplyr package
data.class$lenCat <- as.factor(as.integer(data.class$length %/% 1000))
data.class.byLen <- data.class %>% dplyr::group_by(lenCat, structural_category) %>% dplyr::summarise(count=dplyr::n() ) %>% mutate(perc=count/sum(count))
data.class.byLen$structural_category <- factor(data.class.byLen$structural_category, levels=(xaxislabelsF1), order=TRUE)

p.classByLen.a <- ggplot(data.class.byLen, aes(x=lenCat, y=count, fill=factor(structural_category))) +
  geom_bar(stat='identity', color="black", size=0.3, width=0.85) +
  scale_fill_manual(values = cat.palette, guide='none', name="Structural Category") +
  mytheme+
  theme(legend.position="right") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  scale_y_continuous(expand=c(0,0))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
  labs(x="Transcript length, kb", y="Counts", title="Structural Categories by Transcript Length")


p.classByLen.b <- ggplot(data.class.byLen, aes(x=lenCat, y=perc*100, fill=factor(structural_category))) +
  geom_bar(stat='identity', color ="black", size=0.3, width=0.85) +
  scale_fill_manual(values = cat.palette, guide='none', name="Structural Category") +
  mytheme+
  theme(legend.position="right")  +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  scale_y_continuous(expand=c(0,0))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
  labs(x="Transcript length, kb", y="%", title="Structural Categories by Transcript Length\n\n\n")



##**** PLOT 8: Expression, if isoform expression provided (iso_exp is in TPM)
if (!all(is.na(data.class$iso_exp))){
  p8 <- ggplot(data=data.class, aes(x=structural_category, y=log2(iso_exp+1), fill=structural_category)) +
    geom_boxplot(color="black", size=0.3,  outlier.size = 0.2) +
    scale_x_discrete(drop=FALSE) +
    ylab("log2(TPM+1)") +
    scale_fill_manual(values = cat.palette) +
    guides(fill="none") +
    mytheme  + theme(axis.text.x = element_text(angle = 45)) +
    theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Transcript Expression by Structural Category\n\n" )
}
###Expression, if isoform expression provided (iso_exp is in TPM) by subcategory
if (!all(is.na(data.FSMISM$iso_exp))){
  p8.s1 <- ggplot(data=data.FSMISM, aes(x=subcategory, y=log2(iso_exp+1), fill=subcategory)) +
    geom_boxplot(color="black", size=0.3,  outlier.size = 0.2, position="dodge") +
    scale_x_discrete(drop=TRUE) +
    facet_grid(.~ structural_category, scales = "free_x") +
    ylab("log2(TPM+1)") +
    scale_fill_manual(values = subcat.palette, guide="none") +
    mytheme  + theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position="right", legend.title=element_blank()) +
    theme(axis.text.x  = element_text(size=10))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Transcript Expression by Subcategory\n\n" )
}

if (!all(is.na(data.NICNNC$iso_exp))){
  p8.s2 <- ggplot(data=data.NICNNC, aes(x=subcategory, y=log2(iso_exp+1), fill=subcategory)) +
    geom_boxplot(color="black", size=0.3,  outlier.size = 0.2, position="dodge") +
    scale_x_discrete(drop=TRUE) +
    facet_grid(.~ structural_category, scales = "free_x") +
    ylab("log2(TPM+1)") +
    scale_fill_manual(values = subcat.palette, guide="none") +
    mytheme  + theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position="right", legend.title=element_blank()) +
    theme(axis.text.x  = element_text(size=10))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Transcript Expression by Subcategory\n\n" )
}

if (!all(is.na(data.other$iso_exp))){
  p8.s3 <- ggplot(data=data.other, aes(x=subcategory, y=log2(iso_exp+1), fill=subcategory)) +
    geom_boxplot(color="black", size=0.3,  outlier.size = 0.2, position="dodge") +
    scale_x_discrete(drop=TRUE) +
    facet_grid(.~ structural_category, scales = "free_x") +
    ylab("log2(TPM+1)") +
    scale_fill_manual(values = subcat.palette, guide="none") +
    mytheme  + theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position="right", legend.title=element_blank()) +
    theme(axis.text.x  = element_text(size=10))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Transcript Expression by Subcategory\n\n" )
}



# PLOT 9: FL number, if FL count provided
# convert FL count to TPM
if (!all(is.na(data.class$FL))){
  p9 <- ggplot(data=data.class, aes(x=structural_category, y=log2(FL_TPM+1), fill=structural_category)) +
    geom_boxplot(color="black", size=0.3, outlier.size=0.2) +
    ylab("log2(FL_TPM+1)") +
    scale_x_discrete(drop=FALSE) +
    scale_fill_manual(values = cat.palette) +
    guides(fill="none") +
    mytheme +
    theme(axis.text.x = element_text(angle = 45)) +
    theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Long Reads Count by Structural Category\n\n" )
}

if (!all(is.na(data.FSMISM$FL))){
  p9.s1 <- ggplot(data=data.FSMISM, aes(x=subcategory, y=log2(FL_TPM+1), fill=subcategory)) +
    geom_boxplot(color="black", size=0.3, outlier.size=0.1) +
    facet_grid(.~ structural_category, scales = "free_x") +
    ylab("log2(FL_TPM+1)") +
    scale_x_discrete(drop=TRUE) +
    scale_fill_manual(values = subcat.palette, guide="none") +
    mytheme +
    theme(legend.position="right", legend.title=element_blank()) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(axis.text.x  = element_text(size=10))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Long Reads Count by Subcategory\n\n" )
}

if (!all(is.na(data.NICNNC$FL))){
  p9.s2 <- ggplot(data=data.NICNNC, aes(x=subcategory, y=log2(FL_TPM+1), fill=subcategory)) +
    geom_boxplot(color="black", size=0.3, outlier.size=0.1) +
    facet_grid(.~ structural_category, scales = "free_x") +
    ylab("log2(FL_TPM+1)") +
    scale_x_discrete(drop=TRUE) +
    scale_fill_manual(values = subcat.palette, guide="none") +
    mytheme +
    theme(legend.position="right", legend.title=element_blank()) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(axis.text.x  = element_text(size=10))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Long Reads Count by Subcategory\n\n" )
}

if (!all(is.na(data.other$FL))){
  p9.s3 <- ggplot(data=data.other, aes(x=subcategory, y=log2(FL_TPM+1), fill=subcategory)) +
    geom_boxplot(color="black", size=0.3, outlier.size=0.1) +
    facet_grid(.~ structural_category, scales = "free_x") +
    ylab("log2(FL_TPM+1)") +
    scale_x_discrete(drop=TRUE) +
    scale_fill_manual(values = subcat.palette, guide="none") +
    mytheme +
    theme(legend.position="right", legend.title=element_blank()) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(axis.text.x  = element_text(size=10))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Long Reads Count by Subcategory\n\n" )
}



# PLOT 10: Gene Expression, if expresion provided
if (!all(is.na(data.class$iso_exp))){
  p10 <- ggplot(data=isoPerGene, aes(x=novelGene, y=log2(geneExp+1), fill=novelGene)) +
    geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
    scale_x_discrete(drop=FALSE) +
    xlab("Structural classification") +
    ylab("log2(Gene_TPM+1)") +
    scale_fill_manual(values = myPalette[c(3:4)]) +
    guides(fill="none") +
    mytheme +
    theme(axis.title.x=element_blank()) +
    ggtitle("Annotated vs Novel Gene Expression\n\n" )
}


# PLOT 11: Gene FL number, if FL count provided


if (!all(is.na(data.class$FL))){
  FL_gene <- aggregate(as.integer(data.class$FL), by = list("associatedGene" = data.class$associated_gene), sum)
  colnames(FL_gene)[ncol(FL_gene)] <- "FL_gene"
  isoPerGene <- merge(isoPerGene, FL_gene, by="associatedGene")
  total_fl <- sum(data.class$FL, na.rm=T)
  isoPerGene$FL_gene_TPM <- isoPerGene$FL_gene*(10**6)/total_fl
  
  p11 <- ggplot(data=isoPerGene, aes(x=novelGene, y=log2(FL_gene_TPM+1), fill=novelGene)) +
    geom_boxplot(color="black", size=0.3,outlier.size = 0.2) +
    scale_x_discrete(drop=FALSE) +
    ylab("log2(FL_TPM+1)") +
    scale_fill_manual(values = myPalette[c(3:4)]) +
    guides(fill="none") +
    mytheme +
    theme(axis.title.x=element_blank()) +
    ggtitle("Number of FL reads per Gene by Type of Gene Annotation\n\n" )
  
}



# PLOT 12: NNC expression genes vs not NNC expression genes
# NNC expression genes vs not NNC expression genes

if (!all(is.na(data.class$gene_exp))){
  if (nrow(data.class[data.class$structural_category=="NNC",])!=0){
    
    NNC_genes <- unique(data.class[data.class$structural_category=="NNC","associated_gene"])
    notNNC_genes <- unique(data.class[!data.class$associated_gene%in%NNC_genes,"associated_gene"])
    isoPerGene[isoPerGene$associatedGene %in% notNNC_genes, "NNC_class"] <- "Genes without\n NNC isoforms"
    isoPerGene[isoPerGene$associatedGene %in% NNC_genes, "NNC_class"] <- "Genes with\n NNC isoforms"
    
    isoPerGene$NNC_class <- factor(isoPerGene$NNC_class, levels=c("Genes with\n NNC isoforms","Genes without\n NNC isoforms"),
                                   labels=c("Genes with\n NNC isoforms","Genes without\n NNC isoforms"), order=T)
    
    p12 <- ggplot(data=isoPerGene[!is.na(isoPerGene$NNC_class),], aes(x=NNC_class, y=log2(geneExp+1), fill=NNC_class)) +
      geom_boxplot(color="black", size=0.3, outlier.size=0.2) +
      xlab("") +  
      ylab("log2(Gene_TPM+1)") +
      scale_x_discrete(drop=FALSE) +
      scale_fill_manual(values = c(myPalette[4],"grey38")) +
      guides(fill="none") +
      mytheme +
      theme(axis.title.x=element_blank()) + 
      ggtitle("Gene Expression of NNC And Not NNC Containing Genes\n\n" )
  }
}



if (!all(is.na(data.class$gene_exp))){
  if (nrow(data.class[data.class$structural_category=="NNC",])!=0 & nrow(data.class[data.class$structural_category=="FSM",])!=0 ){
    
    FSM_just_genes = unique(data.class[data.class$FSM_class=="A" & data.class$structural_category=="FSM","associated_gene"])
    NNC_just_genes = unique(data.class[data.class$FSM_class=="A" & data.class$structural_category=="NNC","associated_gene"])
    FSMandNNCgenes = unique(data.class[data.class$FSM_class=="C" & data.class$structural_category=="NNC","associated_gene"])
    isoPerGene[isoPerGene$associatedGene %in% FSMandNNCgenes, "FSM_NNC_class"] <- "Genes expressing\nboth NNC and\n FSM isoforms"
    isoPerGene[isoPerGene$associatedGene %in% NNC_just_genes, "FSM_NNC_class"] <- "Genes expressing\n only NNC isoforms"
    isoPerGene[isoPerGene$associatedGene %in% FSM_just_genes, "FSM_NNC_class"] <- "Genes expressing\n only FSM isoforms"
    data.class[data.class$associated_gene %in% FSMandNNCgenes, "class"] <- "Genes expressing\nboth NNC and\n FSM isoforms"
    data.class[data.class$associated_gene %in% NNC_just_genes, "class"] <- "Genes expressing\n only NNC isoforms"
    data.class[data.class$associated_gene %in% FSM_just_genes, "class"] <- "Genes expressing\n only FSM isoforms"
    
    isoPerGene$FSM_NNC_class = factor(isoPerGene$FSM_NNC_class, levels=c("Genes expressing\nboth NNC and\n FSM isoforms","Genes expressing\n only NNC isoforms","Genes expressing\n only FSM isoforms"),
                                      labels=c("Genes expressing\nboth NNC and\n FSM isoforms","Genes expressing\n only NNC isoforms","Genes expressing\n only FSM isoforms"), order=T)
    
    p13 <- ggplot(data=isoPerGene[!is.na(isoPerGene$FSM_NNC_class),], aes(x=FSM_NNC_class, y=log2(geneExp+1), fill=FSM_NNC_class)) +
      geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
      ylab("log2( # Short reads per gene + 1)") +
      theme(axis.title.x=element_blank()) +
      #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
      scale_fill_manual(values = c("grey38",myPalette[[4]],myPalette[[1]])) +
      guides(fill="none") +
      mytheme +
      theme(axis.title.x=element_blank()) +
      ggtitle("Gene Expression Level in NNC/FSM Containing Genes\n\n" ) +
      scale_x_discrete(breaks=c("Genes expressing\nboth NNC and\n FSM isoforms",
                                "Genes expressing\n only FSM isoforms",
                                "Genes expressing\n only NNC isoforms"),
                       labels=c("NNC/FSM genes",
                                "FSM genes",
                                "NNC genes"), drop=FALSE) 
    
    
    p13.c <- ggplot(data=data.class[!is.na(data.class$class),], aes(x=class, y=log2(iso_exp+1), fill=structural_category)) +
      geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
      ylab("log2( # Short reads per transcript + 1)") +
      theme(axis.title.x=element_blank()) +
      scale_fill_manual(values = myPalette) +
      guides(fill="none") +
      mytheme +
      theme(axis.title.x=element_blank()) +
      ggtitle("Transcript Expression Level in NNC/FSM Containing Genes\n\n" ) +
      scale_x_discrete(breaks=c("Genes expressing\nboth NNC and\n FSM isoforms",
                                "Genes expressing\n only FSM isoforms",
                                "Genes expressing\n only NNC isoforms"),
                       labels=c("NNC/FSM genes",
                                "FSM genes",
                                "NNC genes"), drop=F) 
    
  }
}

if (nrow(data.FSM) > 0) {

  diff_max <- max(max(abs(data.FSM$diff_to_TSS)), max(abs(data.FSM$diff_to_TTS)));
  diff_breaks <- c(-(diff_max+1), seq(-1000, 1000, by = 100), diff_max+1);
  breaks_labels <- c("Larger than -1 kb", "-1 to -0.9 kb", "-0.9 to -0.8 kb", "-0.8 to -0.7 kb", "-0.7 to -0.6 kb", "-0.6 to -0.5 kb",
                     "-0.5 to -0.4 kb", "-0.4 to -0.3 kb", "-0.3 to -0.2 kb", "-0.2 to -0.1 kb", "-0.1 to 0 kb", "0 to 0.1 kb", "0.1 to 0.2 kb",
                     "0.2 to 0.3 kb", "0.3 to 0.4 kb", "0.4 to 0.5 kb", "0.5 to 0.6 kb", "0.6 to 0.7 kb", "0.7 to 0.8 kb", "0.8 to 0.9 kb", "0.9 to 1 kb",
                     "Larger than 1 kb")
  data.FSM$diffTTSCat = cut(-(data.FSM$diff_to_TTS), breaks = diff_breaks);
  data.FSM$diffTSSCat = cut(-(data.FSM$diff_to_TSS), breaks = diff_breaks);
  p21.FSM.list = list()
  p21.FSM.list.a = list()

  # plot histogram of distance to polyA site, Y-axis absolute count
  
  if (!all(is.na(data.FSM$polyA_motif))){
    alpha21='!is.na(polyA_motif)'
    alpha_labs21="polyA motif found"
  }else{
    alpha21=NULL
    alpha_labs21=NULL
  }
  
  max_height <- max(table(data.FSM$diffTTSCat));
  max_height <- (max_height %/% 10+1) * 10;
  p21.a <- ggplot(data=data.FSM, aes(x=diffTTSCat)) +
    geom_bar(fill=myPalette[4], color="black", size=0.3, aes(alpha=eval(parse(text = alpha21)))) +
    scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
    mytheme + labs(alpha = alpha_labs21) +
    scale_x_discrete(drop=F, labels=breaks_labels) +
    ylab("Transcripts, count")+
    xlab("Distance to annotated Transcription Termination Site (TTS), bp")+
    labs(     title="Distance to annotated Transcription Termination Site (TTS)\n\nFSM",
              subtitle="Negative values indicate upstream of annotated termination site\n\n") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
  
  
  p21.b <- ggplot(data=data.FSM, aes(x=diffTTSCat)) +
    geom_bar(aes(y = (..count..)/sum(..count..) , alpha= eval(parse(text = alpha21))), fill=myPalette[4], color="black", size=0.3)+
    scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                       labels=c("0","25","50","75","100"), expand = c(0,0)) +
    scale_x_discrete(drop=F, labels=breaks_labels) +
    mytheme + labs(alpha = alpha_labs21) +
    ylab("Transcripts, %")+
    xlab("Distance to annotated Transcription Termination Site (TTS), bp")+
    labs(     title="Distance to annotated Transcription Termination Site (TTS)\n\nFSM",
              subtitle="Negative values indicate upstream of annotated termination site\n\n") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
  
  p21.stitles.FSM<-list("Distance to annotated Transcription Termination Site (TTS)\n\nFSM Alternative 3'End",
                        "Distance to annotated Transcription Termination Site (TTS)\n\nFSM Alternative 3'5'End",
                        "Distance to annotated Transcription Termination Site (TTS)\n\nFSM Alternative 5'End",
                        "Distance to annotated Transcription Termination Site (TTS)\n\nFSM Reference Match")
  for(i in 1:length(subcategories.FSM)){
    c<-data.frame(subcategories.FSM[i])
    if (!(dim(c))[1]==0 & !all(is.na(c$polyA_motif))){
      diff_max <- max(max(abs(c$diff_to_TSS)), max(abs(c$diff_to_TTS)));
      diff_breaks <- c(-(diff_max+1),seq(-1000, 1000, by = 100),(diff_max+1));
      c$diffTTSCat = cut(-(c$diff_to_TTS), breaks = diff_breaks);
      max_height <- max(table(c$diffTTSCat));
      max_height <- (max_height %/% 10+1) * 10;
      p21.s <- ggplot(data=c, aes(x=diffTTSCat)) +
        geom_bar(fill=myPalette[4], color="black", size=0.3, aes( alpha=eval(parse(text = alpha21)))) +
        scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
        mytheme + labs(alpha = alpha_labs21) +
        scale_x_discrete(drop=F, labels=breaks_labels) +
        ylab("Transcripts, count")+
        xlab("Distance to Annotated Transcription Termination Site (TTS), bp")+
        labs(     title=p21.stitles.FSM[i],
                  subtitle="Negative values indicate upstream of annotated termination site\n\n") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        theme(legend.justification=c(1,1), legend.position=c(1,1))
      p21.s.a <- ggplot(data=c, aes(x=diffTTSCat)) +
        geom_bar(aes(alpha=eval(parse(text = alpha21)), y = (..count..)/sum(..count..)), fill=myPalette[4], color="black", size=0.3)+
        scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                            labels=c("0","25","50","75","100"), expand=c(0,0)) +
        scale_x_discrete(drop=F, labels=breaks_labels) +
        mytheme + labs(alpha = alpha_labs21) +
        ylab("Transcripts, %")+
        xlab("Distance to Annotated Transcription Termination Site (TTS), bp")+
        labs(     title=p21.stitles.FSM[i],
                  subtitle="Negative values indicate upstream of annotated termination site\n\n") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        theme(legend.justification=c(1,1), legend.position=c(1,1))
        
      p21.FSM.list[[i]] = p21.s
      p21.FSM.list.a[[i]] = p21.s.a
    }
    }
  
  # plot histogram of distance to start site, Y-axis absolute count
  max_height <- max(table(data.FSM$diffTSSCat));
  max_height <- (max_height %/% 10+1) * 10;
  if (!all(is.na(data.FSM$within_CAGE_peak))){
    alpha22tss='within_CAGE_peak'
    alpha_labs22tss="TSS within a CAGE peak"
  }else{
    alpha22tss=NULL
    alpha_labs22tss=NULL
  }
  
  p22.a <- ggplot(data=data.FSM, aes(x=diffTSSCat)) +
    geom_bar(fill=myPalette[6], color="black", size=0.3 , aes(alpha=eval(parse(text = alpha22tss))))+
    scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
    scale_x_discrete(drop=F, labels=breaks_labels) +
    mytheme + labs(alpha=alpha_labs22tss) +
    ylab("Transcripts, count") +
    xlab("Distance to Annotated Transcription Start Site (TSS), bp") +
    labs(title="Distance to Annotated Transcription Start Site for FSM\n\n",
         subtitle="Negative values indicate downstream of annotated TSS\n\n") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
  
  p22.b <- ggplot(data=data.FSM, aes(x=diffTSSCat)) +
    geom_bar(aes(y = (..count..)/sum(..count..) , alpha=eval(parse(text = alpha22tss))), fill=myPalette[6], color="black", size=0.3)+
    scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                       labels=c("0","25","50","75","100"), expand = c(0,0)) +
    scale_x_discrete(drop=F, labels=breaks_labels) +
    mytheme + labs(alpha=alpha_labs22tss) +
    ylab("Transcripts, %")+
    xlab("Distance to Annotated Transcription Start Site (TSS), bp")+
    labs(title="Distance to Annotated Transcription Start Site for FSM\n\n",
         subtitle="Negative values indicate downstream of annotated TSS\n\n") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
  
    #FSM_TSS
  p22.FSM.list = list()
  p22.FSM.list.a = list()
  p22.stitles.FSM<-list("Distance to Annotated Transcription Start Site\n\nFSM Alternative 3' End",
                        "Distance to Annotated Transcription Start Site\n\nFSM Alternative 3'5' End",
                        "Distance to Annotated Transcription Start Site\n\nFSM Alternative 5' End",
                        "Distance to Annotated Transcription Start Site\n\nFSM Reference Match")
    
  for(i in 1:length(subcategories.FSM)){
    c<-data.frame(subcategories.FSM[i])
    if (!(dim(c))[1]==0 & !all(is.na(c$within_CAGE_peak))){
      diff_max <- max(max(abs(c$diff_to_TSS)), max(abs(c$diff_to_TTS)));
      diff_breaks <- c(-(diff_max+1),seq(-1000, 1000, by = 100), (diff_max+1));
      c$diffTTSCat = cut(-(c$diff_to_TTS), breaks = diff_breaks);
      c$diffTSSCat = cut(-(c$diff_to_TSS), breaks = diff_breaks);
      max_height <- max(table(c$diffTSSCat));
      max_height <- (max_height %/% 10+1) * 10;
      p22.s <- ggplot(data=c, aes(x=diffTSSCat)) +
        geom_bar(fill=myPalette[6], color="black", size=0.3, aes(alpha=eval(parse(text = alpha22tss)))) +
        scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
        scale_x_discrete(drop=F, labels=breaks_labels) +
        mytheme + labs(alpha=alpha_labs22tss) +
        ylab("Transcripts, count")+
        xlab("Distance to annotated Transcription Start Site (TSS), bp")+
        labs(     title=p22.stitles.FSM[i],
                  subtitle="Negative values indicate downstream of annotated TSS\n\n") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        theme(legend.justification=c(1,1), legend.position=c(1,1))
      p22.s.a <- ggplot(data=c, aes(x=diffTSSCat)) +
        geom_bar(aes(alpha=eval(parse(text = alpha22tss)), y = (..count..)/sum(..count..)), fill=myPalette[6], color="black", size=0.3)+
        scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                            labels=c("0","25","50","75","100"), expand = c(0,0)) +
        scale_x_discrete(drop=F, labels=breaks_labels) +
        mytheme + labs(alpha=alpha_labs22tss) +
        ylab("Transcripts, %")+
        xlab("Distance to annotated Transcription Start Site (TSS), bp")+
        labs(     title=p22.stitles.FSM[i],
                  subtitle="Negative values indicate downstream of annotated TSS\n\n") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        theme(legend.justification=c(1,1), legend.position=c(1,1))
      p22.FSM.list[[i]] = p22.s
      p22.FSM.list.a[[i]] = p22.s.a
    }
  }
}


  
if (nrow(data.ISM) > 0) {
  
  diff_max <- max(max(abs(data.ISM$diff_to_TSS)), max(abs(data.ISM$diff_to_TTS)));
  diff_breaks <- c(-(diff_max+1), seq(-5000, 5000, by = 500), diff_max+1);
  breaks_labels <- c("Larger than -5 kb", "-5 to -4.5 kb", "-4.5 to -4 kb", "-4 to -3.5 kb", "-3.5 to -3 kb", "-3 to -2.5 kb",
                     "-2.5 to -2 kb", "-2 to -1.5 kb", "-1.5 to -1 kb", "-1 to -0.5 kb", "-0.5 to 0 kb", "0 to 0.5 kb", "0.5 to 1 kb",
                     "1 to 1.5 kb", "1.5 to 2 kb", "2 to 2.5 kb", "2.5 to 3 kb", "3 to 3.5 kb", "3.5 to 4 kb", "4 to 4.5 kb", "4.5 to 5 kb", "Larger than 5 kb")
  
  data.ISM$diffTTSCat = cut(-(data.ISM$diff_to_TTS), breaks = diff_breaks);
  data.ISM$diffTSSCat = cut(-(data.ISM$diff_to_TSS), breaks = diff_breaks);
  
  max_height <- max(table(data.ISM$diffTTSCat));
  max_height <- (max_height %/% 10+1) * 10;
  
  p21.dist3.ISM.a <- ggplot(data=data.ISM, aes(x=diffTTSCat)) +
    geom_bar(fill=myPalette[4], color="black", size=0.3,  aes( alpha= !is.na(polyA_motif))) +
    scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
    mytheme + labs(alpha = "polyA motif found") +
    scale_x_discrete(drop=F, labels=breaks_labels) +
    ylab("Transcripts, count")+
    xlab("Distance to annotated polyadenylation site, bp")+
    labs(     title="Distance to Annotated Polyadenylation Site for ISM\n\n",
              subtitle="Negative values indicate upstream of annotated polyA site\n\n") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
  
  # plot histogram of distance to polyA site, Y-axis percentages
  p21.dist3.ISM.b <- ggplot(data=data.ISM, aes(x=diffTTSCat)) +
    geom_bar(aes(y = (..count..)/sum(..count..), alpha= !is.na(polyA_motif)), fill=myPalette[4], color="black", size=0.3)+
    scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                       labels=c("0","25","50","75","100"), expand = expansion(mult=c(0,0.1))) +
    scale_x_discrete(drop=F, labels=breaks_labels) +
    mytheme + labs(alpha = "polyA motif found") +
    ylab("Transcripts, %")+
    xlab("Distance to annotated polyadenylation site, bp")+
    labs(     title="Distance to Annotated Polyadenylation Site for ISM\n\n",
              subtitle="Negative values indicate upstream of annotated polyA site\n\n") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
  
  # plot histogram of distance to start site, Y-axis absolute count
  p22.dist5.ISM.a <- ggplot(data=data.ISM, aes(x=diffTSSCat)) +
    geom_bar(fill=myPalette[6], color="black", size=0.3, aes(alpha=within_CAGE_peak))+
    scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
    scale_x_discrete(drop=F, labels=breaks_labels) +
    mytheme + labs(alpha="TSS within a CAGE peak") +
    ylab("Transcripts, count")+
    xlab("Distance to annotated transcription start site, bp")+
    labs(     title="Distance to Annotated Transcription Start Site for ISM\n\n",
              subtitle="Negative values indicate downstream of annotated TSS\n\n") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
  
  # plot histogram of distance to start site, Y-axis absolute count
  p22.dist5.ISM.b <- ggplot(data=data.ISM, aes(x=diffTSSCat)) +
    geom_bar(aes(y = (..count..)/sum(..count..), alpha=within_CAGE_peak), fill=myPalette[6], color="black", size=0.3)+
    scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                       labels=c("0","25","50","75","100"), expand = c(0,0)) +
    scale_x_discrete(drop=F, labels=breaks_labels) +
    mytheme + labs(alpha="TSS within a CAGE peak") +
    ylab("Transcripts, %")+
    xlab("Distance to annotated transcription start site, bp")+
    labs(title="Distance to Annotated Transcription Start Site for ISM\n\n",
         subtitle="Negative values indicate downstream of annotated TSS\n\n") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    theme(legend.justification=c(1,1), legend.position=c(1,1))
  p21.stitles.ISM<-list("Distance to Annotated Polyadenylation Site for ISM\n\n3' Fragment",
                        "Distance to Annotated Polyadenylation Site for ISM\n\nInternal Fragment",
                        "Distance to Annotated Polyadenylation Site for ISM\n\nA5' Fragment",
                        "Distance to Annotated Polyadenylation Site for ISM\n\nIntron Retention")
  p21.ISM.list = list()
  p21.ISM.list.a = list()
  for(i in 1:length(subcategories.ISM)){
    c<-data.frame(subcategories.ISM[i])
    if (!(dim(c))[1]==0 & !all(is.na(c$polyA_motif))){
      diff_max <- max(max(abs(c$diff_to_TSS)), max(abs(c$diff_to_TTS)));
      c$diffTTSCat = cut(-(c$diff_to_TTS), breaks = diff_breaks);
      max_height <-  max(table(c$diffTTSCat));
      max_height <- (max_height %/% 10+1) * 10;
      p21.s <- ggplot(data=c, aes(x=diffTTSCat)) +
        geom_bar(fill=myPalette[4], color="black", size=0.3, aes( alpha= !is.na(polyA_motif))) +
        scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
        scale_x_discrete(drop=F, labels=breaks_labels) +
        mytheme + labs(alpha = "polyA motif found") +
        ylab("Transcripts, count")+
        xlab("Distance to annotated polyadenylation site, bp")+
        labs(     title=p21.stitles.ISM[i],
                  subtitle="Negative values indicate upstream of annotated polyA site\n\n") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        theme(legend.justification=c(1,1), legend.position=c(1,1))
      p21.s.a <- ggplot(data=c, aes(x=diffTTSCat)) +
        geom_bar(aes(alpha= !is.na(polyA_motif), y = (..count..)/sum(..count..)), fill=myPalette[4], color="black", size=0.3)+
        scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                           labels=c("0","25","50","75","100"), expand = c(0,0)) +
        scale_x_discrete(drop=F, labels=breaks_labels) +
        mytheme + labs(alpha = "polyA motif found") +
        ylab("Transcripts, %")+
        xlab("Distance to annotated polyadenylation site, bp")+
        labs(     title=p21.stitles.ISM[i],
                  subtitle="Negative values indicate upstream of annotated polyA site\n\n") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        theme(legend.justification=c(1,1), legend.position=c(1,1))
      
      p21.ISM.list[[i]] = p21.s
      p21.ISM.list.a[[i]] = p21.s.a
    }
  }
  
  #FSM_TSS
  p22.stitles.ISM<-list("Distance to Annotated Transcription Start Site for ISM\n\n3' Fragment",
                        "Distance to Annotated Transcription Start Site for ISM\n\nInternal Fragment",
                        "Distance to Annotated Transcription Start Site for ISM\n\nA5' Fragment",
                        "Distance to Annotated Transcription Start Site for ISM\n\nIntron Retention")
  p22.ISM.list = list()
  p22.ISM.list.a = list()
  for(i in 1:length(subcategories.ISM)){
    c<-data.frame(subcategories.ISM[i])
    if (!(dim(c))[1]==0 & !all(is.na(c$within_CAGE_peak))){
      c$diffTTSCat = cut(-(c$diff_to_TTS), breaks = diff_breaks);
      c$diffTSSCat = cut(-(c$diff_to_TSS), breaks = diff_breaks);
      max_height <- max(max(table(c$diffTSSCat)), max(table(c$diffTTSCat)));
      max_height <- (max_height %/% 10+1) * 10;
      p22.s <- ggplot(data=c, aes(x=diffTSSCat)) +
        geom_bar(fill=myPalette[6], color="black", size=0.3, aes( alpha= within_CAGE_peak)) +
        scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
        scale_x_discrete(drop=F, labels=breaks_labels) +
        mytheme + labs(alpha = "TSS within a CAGE peak") +
        ylab("Transcripts, count")+
        xlab("Distance to annotated transcription start site, bp")+
        labs(     title=p22.stitles.ISM[i],
                  subtitle="Negative values indicate downstream of annotated TSS\n\n") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        theme(legend.justification=c(1,1), legend.position=c(1,1))
      p22.s.a <- ggplot(data=c, aes(x=diffTSSCat)) +
        geom_bar(aes( alpha= within_CAGE_peak, y = (..count..)/sum(..count..)), fill=myPalette[6], color="black", size=0.3)+
        scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                           labels=c("0","25","50","75","100"), expand = c(0,0)) +
        scale_x_discrete(drop=F, labels=breaks_labels) +
        mytheme + labs(alpha = "TSS within a CAGE peak") +
        ylab("Transcripts, %")+
        xlab("Distance to annotated transcription start site, bp")+
        labs(     title=p22.stitles.ISM[i],
                  subtitle="Negative values indicate downstream of annotated TSS\n\n") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        theme(legend.justification=c(1,1), legend.position=c(1,1))
      p22.ISM.list[[i]] = p22.s
      p22.ISM.list.a[[i]] = p22.s.a
    }
  }
}



# PLOT 23: Junction categories


if (nrow(data.junction) > 0){
  data.junction$junctionLabel = with(data.junction, paste(chrom, strand,genomic_start_coord, genomic_end_coord, sep="_"))
  
  data.junction$canonical_known = with(data.junction, paste(junction_category,canonical,"SJ", sep="_"))
  data.junction$canonical_known=as.factor(data.junction$canonical_known)
  data.junction$canonical_known = factor(data.junction$canonical_known, levels=c("known_canonical_SJ", "known_non_canonical_SJ", "novel_canonical_SJ", "novel_non_canonical_SJ"),
                                         labels=c("Known\ncanonical ", "Known\nNon-canonical ", "Novel\ncanonical ", "Novel\nNon-canonical "), order=T) 
  data.junction$structural_category = data.class[data.junction$isoform,"structural_category"]
  ##    data.junction$TSSrange =cut(data.junction$transcript_coord, breaks = c(0, 40, 80, 120, 160, 200, 10000000), labels = c("0-40", "41-80", "81-120", "121-160", "161-200",">200"))
  
  p23.a <- ggplot(data.junction, aes(x=structural_category)) +
    geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=SJ_type), color="black",  size=0.3, width = 0.7) +
    scale_y_continuous(breaks=c(0.0,0.25,0.5,0.75,1),
                       labels=c("0","25","50","75","100"), expand = c(0,0)) +
    scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
    ylab("Splice junctions, %") +
    mytheme +
    guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.3))+
    theme(legend.position="bottom", legend.title=element_blank())  +
    theme(axis.text.x = element_text(angle = 45)) +
    theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Distribution of Splice Junctions by Structural Classification\n\n\n")
  
  t <- subset(data.class, exons > 1)  # select only multi-exon isoforms
  
  p23.b <- ggplot(data=t, aes(x=structural_category)) +
    geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=all_canonical), color="black", size=0.3, width = 0.7) +
    scale_y_continuous(breaks=c(0.0,0.25,0.5,0.75,1),
                       labels=c("0","25","50","75","100"), expand = c(0,0)) +
    scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
    xlab("") +
    ylab("Transcripts, %") +
    mytheme +
    guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.3))+
    theme(legend.position="bottom", legend.title=element_blank())  +
    theme(axis.text.x = element_text(angle = 45)) +
    theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
    theme(axis.title.x=element_blank()) +
    ggtitle("Distribution of Transcripts by Splice Junctions\n\n\n")
}

# PLOT 29: RT-switching


if (sum(data.junction$RTS_junction=='TRUE') > 0) {
  
  a <- data.frame(table(data.junction$SJ_type));
  b <- data.frame(table(subset(data.junction, RTS_junction=="TRUE")$SJ_type));
  
  df.RTS <- merge(a, b, by="Var1");
  df.RTS$perc <- df.RTS$Freq.y/df.RTS$Freq.x *100
  df.RTS[is.na(df.RTS$perc),"perc"] <- 0
  
  maxH <- min(100, (max(df.RTS$perc) %/% 5) * 5 + 5);
  
  p29.a <- ggplot(data=df.RTS, aes(x=Var1, y=perc, fill=Var1)) +
    geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
    geom_text(label=paste(round(df.RTS$perc, 2),"%",sep=''), nudge_y=0.3) +
    scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=F) +
    labs(x="", y="RT-switching junctions, %") +
    ggtitle("RT-Switching All Junctions\n\n" ) +
    mytheme +
    guides(fill="none") +
    scale_y_continuous(expand = c(0,0), limits = c(0,maxH)) +
    theme(axis.text.x = element_text(size=11))
  
  c <- data.frame(table(uniqJuncRTS$SJ_type));
  d <- data.frame(table(subset(uniqJuncRTS, RTS_junction=='TRUE')$SJ_type));
  
  
  df.uniqRTS <- merge(c, d, by="Var1");
  df.uniqRTS$perc <- df.uniqRTS$Freq.y/df.uniqRTS$Freq.x *100
  df.uniqRTS[is.na(df.uniqRTS$perc),"perc"] <- 0
  
  maxH <- min(100, (max(df.uniqRTS$perc) %/% 5) * 5 + 5);
  
  p29.b <- ggplot(data=df.uniqRTS, aes(x=Var1, y=perc, fill=Var1)) +
    geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
    geom_text(label=paste(round(df.uniqRTS$perc, 2),"%",sep=''), nudge_y=0.3) +
    scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=F) +
    labs(x="", y="RT-switching junctions, %") +
    ggtitle("Unique Junctions RT-switching\n\n" ) +
    mytheme +
    guides(fill="none") +
    scale_y_continuous(expand = c(0,0), limits = c(0,maxH)) +
    theme(axis.text.x = element_text(size=11))
  
}

# PLOT pn4-5: Splice Junction Coverage (if coverage provided)

if (!all(is.na(data.junction$total_coverage_unique))){
  
  uniqJuncCov <- unique(data.junction[,c("junctionLabel","SJ_type", "total_coverage_unique")])
  
  e <- data.frame(table(uniqJuncCov$SJ_type))
  f <- data.frame(table(uniqJuncCov[which(uniqJuncCov$total_coverage_unique>0),"SJ_type"]))
  
  
  df.juncSupport <- data.frame(type=e$Var1, count=e$Freq-f$Freq, name='Unsupported')
  df.juncSupport <- rbind(df.juncSupport, data.frame(type=f$Var1, count=f$Freq, name='Supported'))
  
  pn4.a <- ggplot(df.juncSupport, aes(x=type, y=count, fill=name)) +
    geom_bar(stat='identity', color="black", size=0.3, width=0.7) +
    scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
    scale_y_continuous( expand = c(0,0)) +
    labs(x='', y='Junctions, count', title='Unique Junctions w/ or w/out Short Read Coverage\n\n\n') +
    mytheme +
    theme(legend.position="bottom", legend.title=element_blank()) +
    guides(fill = guide_legend(title = "") )
  
  
  df.SJcov <- merge(e, f, by="Var1")
  # calculate the percentage of junctions that have zero short read junction coverage
  df.SJcov$perc <- 100-df.SJcov$Freq.y / df.SJcov$Freq.x * 100 ;
  df.SJcov[is.na(df.SJcov$perc), "perc"] <- 0
  maxP=min(100, (max(df.SJcov$perc) %/% 5) * 5) + 5 ;
  pn4.b <- ggplot(df.SJcov, aes(x=Var1,fill=Var1, y=perc)) +
    geom_bar(stat="identity", position = position_dodge(), color="black", size=0.3, width=0.7) +
    geom_text(label=paste(round(df.SJcov$perc,1.05),"%",sep=''), nudge_y=1.5) +
    scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
    scale_y_continuous( expand = c(0,0), limits = c(0,maxP)) +
    labs(x='', y="Junctions, %", title='Unique Junctions w/out Short Read Coverage\n\n\n') +
    mytheme +
    guides(fill="none")
  
}

##### Distances to CAGE peaks by FSM and ISM

if (!all(is.na(data.class$dist_to_CAGE_peak))) {
  diff_max=11000     ##diff_max <- max(abs(data.class$dist_to_cage_peak));
  diff_breaks <- c(-(diff_max+1), seq(-200, 100, by = 20), diff_max+1);
  
  breaks_labels <- c("Larger than -200", "-200 to -180","-180 to -160","-160 to -140","-140 to -120",
                     "-120 to -100", "-100 to -80", "-80 to -60", "-60 to -40", "-40 to -20", "-20 to 0",
                     "0 to 20", "20 to 40", "40 to 60", "60 to 80", "80 to 100", "Larger than 100")
  
  data.FSM$dist_CAGE_Cat = cut(-(data.FSM$dist_to_CAGE_peak), breaks = diff_breaks);
  data.ISM$dist_CAGE_Cat = cut(-(data.ISM$dist_to_CAGE_peak), breaks = diff_breaks);
  data.NNC$dist_CAGE_Cat = cut(-(data.NNC$dist_to_CAGE_peak), breaks = diff_breaks);
  data.NIC$dist_CAGE_Cat = cut(-(data.NIC$dist_to_CAGE_peak), breaks = diff_breaks);
  
  #  max_height <- max(max(table(data.class$dist_cage_Cat)), max(table(data.class$dist_cage_Cat)));
  #  max_height <- (max_height %/% 10+1) * 10;
  
  #  ggplot(data=d.fsm, aes(x=dist_to_cage_peak , fill=structural_category)) +
  #    geom_density( color="black", size=0.3) + xlim(c(-50,50))+
  #   scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
  #    mytheme  + facet_wrap(~structural_category , nrow=2)+
  #   scale_x_discrete(limits = c(-50,50)) +
  #    ylab("Number of transcripts")+
  #    xlab("Distance to CAGE peak, bp")+
  #    labs(     title="Distance to CAGE Peak",
  #            subtitle="Negative values indicate downstream of annotated CAGE peak\n\n") +
  #    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  ## FSM cage hist number of Isoforms  
  cage_hist_FSM=ggplot(data=subset(data.FSM, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat, fill=structural_category)) +
    geom_bar(aes(alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[1]) +
    mytheme  +
    scale_x_discrete(drop=F, labels=breaks_labels) +
    scale_y_continuous( expand = expansion(mult = c(0,0.1)) ) +
    theme(legend.position="bottom") +
    ylab("Number of transcripts")+
    xlab("Distance to CAGE peak, bp")+
    labs(     title="Distance to CAGE Peak of Multi-Exonic FSM ",
              subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
              alpha = "TSS within a CAGE peak") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  cage_hist_FSM_perc=ggplot(data=subset(data.FSM, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
    geom_bar(aes(y = (..count..)/sum(..count..), alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[1]) +
    scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                       labels=c("0","25","50","75","100"), 
                       expand = expansion(mult = c(0,0.1))) +
    mytheme  +
    theme(legend.position="bottom") +
    scale_x_discrete(drop=F, labels=breaks_labels) +
    ylab("Transcripts, %")+
    xlab("Distance to CAGE peak, bp")+
    labs(     title="Distance to CAGE Peak of Multi-Exonic FSM",
              subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
              alpha = "TSS within a CAGE peak") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  cage.titles.FSM<-list("Distance to CAGE Peak of Multi-Exonic FSM\n\nAlternative 3' End",
                        "Distance to CAGE Peak of Multi-Exonic FSM\n\nAlternative 3'5' End",
                        "Distance to CAGE Peak of Multi-Exonic FSM\n\nAlternative 5' End",
                        "Distance to CAGE Peak of Multi-Exonic FSM\n\nReference Match")
  cage.FSM.list = list()
  cage.FSM.list.a = list()
  for(i in 1:length(subcategories.FSM)){
    c<-data.frame(subcategories.FSM[i])
    if (!(dim(c))[1]==0 & !all(is.na(c$dist_to_CAGE_peak))){
      #diff_max=11000
      #diff_breaks <- c(-(diff_max+1), seq(-200, 100, by = 20), diff_max+1);
      #c$dist_to_cage_peak[is.na(c$dist_to_cage_peak)] <- -11000
      c$dist_CAGE_Cat = cut(-(c$dist_to_CAGE_peak), breaks = diff_breaks);
      cage.FSM.s=ggplot(data=subset(c, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
        geom_bar(aes(alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[1]) +
        mytheme  +
        scale_x_discrete(drop=F, labels=breaks_labels) +
        scale_y_continuous( expand = expansion(mult = c(0,0.1)) ) +
        theme(legend.position="bottom") +
        ylab("Number of transcripts")+
        xlab("Distance to CAGE peak, bp")+
        labs(     title=cage.titles.FSM[i],
                  subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                  alpha = "TSS within a CAGE peak") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
      cage.FSM.s.perc=ggplot(data=subset(c, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
        geom_bar(aes(y = (..count..)/sum(..count..), alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[1]) +
        scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                           labels=c("0","25","50","75","100"), 
                           expand = expansion(mult = c(0,0.1))) +
        mytheme  +
        scale_x_discrete(drop=F, labels=breaks_labels) +
        theme(legend.position="bottom") +
        ylab("Transcripts, %")+
        xlab("Distance to CAGE peak, bp")+
        labs(     title=cage.titles.FSM[i],
                  subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                  alpha = "TSS within a CAGE peak") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
      cage.FSM.list[[i]] = cage.FSM.s
      cage.FSM.list.a[[i]] = cage.FSM.s.perc
    }
  }
  
  cage_hist_ISM=ggplot(data=subset(data.ISM, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
    geom_bar(aes(alpha=within_CAGE_peak), color="black", size=0.3, fill=myPalette[2]) +
    #    scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
    mytheme  +
    scale_x_discrete(drop=F, labels=breaks_labels) +
    scale_y_continuous( expand = expansion(mult = c(0,0.1)) ) +
    theme(legend.position="bottom") +
    #    scale_x_discrete(limits = c(-50,50)) +
    ylab("Number of transcripts")+
    xlab("Distance to CAGE peak, bp")+
    labs(     title="Distance to CAGE Peak of Multi-Exonic ISM",
              subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
              alpha = "TSS within a CAGE peak") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  cage_hist_ISM_perc=ggplot(data=subset(data.ISM, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
    geom_bar(aes(y = (..count..)/sum(..count..), alpha=within_CAGE_peak), color="black", size=0.3, fill=myPalette[2]) +
    scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                       labels=c("0","25","50","75","100"), 
                       expand = expansion(mult = c(0,0.1))) +
    scale_x_discrete(drop=F, labels=breaks_labels) +
    mytheme  + theme(legend.position="bottom") +
    #    scale_x_discrete(limits = c(-50,50)) +
    ylab("Transcripts, %")+
    xlab("Distance to CAGE Peak, bp")+
    labs(     title="Distance to CAGE Peak of Multi-Exonic ISM",
              subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
              alpha = "TSS within a CAGE peak") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  cage.titles.ISM<-list("Distance to CAGE Peak of Multi-Exonic ISM\n\n3' Fragment",
                        "Distance to CAGE Peak of Multi-Exonic ISM\n\nInternal Fragment",
                        "Distance to CAGE Peak of Multi-Exonic ISM\n\n5' Fragment",
                        "Distance to CAGE Peak of Multi-Exonic ISM\n\nIntron Retention")
  
  
  cage.ISM.list = list()
  cage.ISM.list.a = list()
  for(i in 1:length(subcategories.ISM)){
    c<-data.frame(subcategories.ISM[i])
    if (!(dim(c))[1]==0 & !all(is.na(c$dist_to_CAGE_peak))){
      #diff_max=11000
      #diff_breaks <- c(-(diff_max+1), seq(-500, 500, by = 20), diff_max+1);
      #c$dist_to_cage_peak[is.na(c$dist_to_cage_peak)] <- -11000
      c$dist_CAGE_Cat = cut(-(c$dist_to_CAGE_peak), breaks = diff_breaks);
      cage.ISM.s=ggplot(data=subset(c, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
        geom_bar(aes(alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[2]) +
        mytheme  +
        scale_x_discrete(drop=F, labels=breaks_labels) +
        scale_y_continuous( expand = expansion(mult = c(0,0.1)) ) +
        theme(legend.position="bottom") +
        ylab("Number of transcripts")+
        xlab("Distance to CAGE peak, bp")+
        labs(     title=cage.titles.ISM[i],
                  subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                  alpha = "TSS within a CAGE peak") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
      cage.ISM.s.perc=ggplot(data=subset(c, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
        geom_bar(aes(y = (..count..)/sum(..count..), alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[2]) +
        scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                           labels=c("0","25","50","75","100"), 
                           expand = expansion(mult = c(0,0.1))) +
        mytheme  +
        scale_x_discrete(drop=F, labels=breaks_labels) +
        theme(legend.position="bottom") +
        ylab("Transcripts, %")+
        xlab("Distance to CAGE peak, bp")+
        labs(     title=cage.titles.ISM[i],
                  subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                  alpha = "TSS within a CAGE peak") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
      cage.ISM.list[[i]] = cage.ISM.s
      cage.ISM.list.a[[i]] = cage.ISM.s.perc
    }
  }
  
  cage_hist_NIC=ggplot(data=subset(data.NIC, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
    geom_bar( aes(alpha=within_CAGE_peak),color="black", size=0.3, fill=myPalette[3]) +
    #    scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
    mytheme  + theme(legend.position="bottom") +
    scale_x_discrete(drop=F, labels=breaks_labels) +
    scale_y_continuous( expand = expansion(mult = c(0,0.1)) ) +
    ylab("Number of transcripts")+
    xlab("Distance to CAGE peak, bp")+
    labs(     title="Distance to CAGE Peak of Multi-Exonic NIC",
              subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
              alpha = "TSS within a CAGE peak") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  cage_hist_NIC_perc=ggplot(data=subset(data.NIC, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
    geom_bar( aes(y = (..count..)/sum(..count..),alpha=within_CAGE_peak),color="black", size=0.3, fill=myPalette[3]) +
    scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                       labels=c("0","25","50","75","100"), 
                       expand = expansion(mult = c(0,0.1))) +
    mytheme  + theme(legend.position="bottom") +
    scale_x_discrete(drop=F, labels=breaks_labels) +
    #    scale_x_discrete(limits = c(-50,50)) +
    ylab("Transcripts, %")+
    xlab("Distance to CAGE peak, bp")+
    labs(     title="Distance to CAGE Peak of Multi-Exonic NIC",
              subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
              alpha = "TSS within a CAGE peak") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  cage.titles.NIC<-list("Distance to CAGE Peak of Multi-Exonic NIC\n\nCombination of Annotated Junctions",
                        "Distance to CAGE Peak of Multi-Exonic NIC\n\nCombination of Annotated Splice Sites",
                        "Distance to CAGE Peak of Multi-Exonic NIC\n\nIntron Retention",
                        "Distance to CAGE Peak of Multi-Exonic NIC\n\nMono-Exon by Intron Retention")
  cage.NIC.list = list()
  cage.NIC.list.a = list()
  for(i in 1:length(subcategories.NIC)){
    c<-data.frame(subcategories.NIC[i])
    if (!(dim(c))[1]==0 & !all(is.na(c$dist_to_CAGE_peak))){
      #diff_max=11000
      #diff_breaks <- c(-(diff_max+1), seq(-500, 500, by = 20), diff_max+1);
      #c$dist_to_cage_peak[is.na(c$dist_to_cage_peak)] <- -11000
      c$dist_CAGE_Cat = cut(-(c$dist_to_CAGE_peak), breaks = diff_breaks);
      cage.NIC.s=ggplot(data=subset(c, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
        geom_bar(aes(alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[3]) +
        mytheme  +
        scale_x_discrete(drop=F, labels=breaks_labels) +
        scale_y_continuous( expand = expansion(mult = c(0,0.1)) ) +
        theme(legend.position="bottom") +
        ylab("Number of transcripts")+
        xlab("Distance to CAGE peak, bp")+
        labs(     title=cage.titles.NIC[i],
                  subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                  alpha = "TSS within a CAGE peak") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
      cage.NIC.s.perc=ggplot(data=subset(c, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
        geom_bar(aes(y = (..count..)/sum(..count..), alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[3]) +
        scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                           labels=c("0","25","50","75","100"), 
                           expand = expansion(mult = c(0,0.1))) +
        mytheme  +
        scale_x_discrete(drop=F, labels=breaks_labels) +
        theme(legend.position="bottom") +
        ylab("Transcripts, %")+
        xlab("Distance to CAGE peak, bp")+
        labs(     title=cage.titles.NIC[i],
                  subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                  alpha = "TSS within a CAGE peak") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
      cage.NIC.list[[i]] = cage.NIC.s
      cage.NIC.list.a[[i]] = cage.NIC.s.perc
    }
  }
  
  cage_hist_NNC=ggplot(data=subset(data.NNC, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
    geom_bar(aes(alpha=within_CAGE_peak), color="black", size=0.3, fill=myPalette[4]) +
    #    scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
    mytheme  + theme(legend.position="bottom") +
    #    scale_x_discrete(limits = c(-50,50)) +
    scale_x_discrete(drop=F, labels=breaks_labels) +
    ylab("Number of transcripts")+
    xlab("Distance to CAGE peak, bp")+
    labs(     title="Distance to CAGE Peak of Multi-Exonic NNC",
              subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
              alpha = "TSS within a CAGE peak") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  cage_hist_NNC_perc=ggplot(data=subset(data.NNC, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
    geom_bar(aes(y = (..count..)/sum(..count..),alpha=within_CAGE_peak), color="black", size=0.3, fill=myPalette[4]) +
    scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                       labels=c("0","25","50","75","100"), 
                       expand = expansion(mult = c(0,0.1))) +
    mytheme + theme(legend.position="bottom") +
    scale_x_discrete(drop=F, labels=breaks_labels) +
    #    scale_x_discrete(limits = c(-50,50)) +
    ylab("Transcripts, %")+
    xlab("Distance to CAGE peak, bp")+
    labs(     title="Distance to CAGE Peak of Multi-Exonic NNC",
              subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
              alpha = "TSS within a CAGE peak") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  cage.titles.NNC<-list("Distance to CAGE Peak of Multi-Exonic NNC\n\nCombination of Annotated Junctions",
                        "Distance to CAGE Peak of Multi-Exonic NNC\n\nCombination of Annotated Splice Sites",
                        "Distance to CAGE Peak of Multi-Exonic NNC\n\nIntron Retention",
                        "Distance to CAGE Peak of Multi-Exonic NNC\n\nMono-Exon by Intron Retention",
                        "Distance to CAGE Peak of Multi-Exonic NNC\n\nAt Least One Annotated Donor/Acceptor")
  
  cage.NNC.list = list()
  cage.NNC.list.a = list()
  for(i in 1:length(subcategories.NNC)){
    c<-data.frame(subcategories.NNC[i])
    if (!(dim(c))[1]==0 & !all(is.na(c$dist_to_CAGE_peak))){
      #diff_max=11000
      #diff_breaks <- c(-(diff_max+1), seq(-500, 500, by = 20), diff_max+1);
      #c$dist_to_cage_peak[is.na(c$dist_to_cage_peak)] <- -11000
      c$dist_CAGE_Cat = cut(-(c$dist_to_CAGE_peak), breaks = diff_breaks);
      cage.NNC.s=ggplot(data=subset(c, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
        geom_bar(aes(alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[4]) +
        mytheme  +
        scale_x_discrete(drop=F, labels=breaks_labels) +
        scale_y_continuous( expand = expansion(mult = c(0,0.1)) ) +
        theme(legend.position="bottom") +
        ylab("Number of transcripts")+
        xlab("Distance to CAGE peak, bp")+
        labs(     title=cage.titles.NNC[i],
                  subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                  alpha = "TSS within a CAGE peak") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
      cage.NNC.s.perc=ggplot(data=subset(c, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
        geom_bar(aes(y = (..count..)/sum(..count..), alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[4]) +
        scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                           labels=c("0","25","50","75","100"), 
                           expand = expansion(mult = c(0,0.1))) +
        mytheme  +
        scale_x_discrete(drop=F, labels=breaks_labels) +
        theme(legend.position="bottom") +
        ylab("Transcripts, %")+
        xlab("Distance to CAGE peak, bp")+
        labs(     title=cage.titles.NNC[i],
                  subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                  alpha = "TSS within a CAGE peak") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
      cage.NNC.list[[i]] = cage.NNC.s
      cage.NNC.list.a[[i]] = cage.NNC.s.perc
    }
  }
  
  

}

#########
if (sum(!is.na(data.class$polyA_dist)) > 10) {
p.polyA_dist <- ggplot(data.class, aes(x=polyA_dist, color=structural_category)) +
    geom_freqpoly(binwidth=1, size=1) +
    scale_color_manual(values = cat.palette)+
    xlab("Distance of polyA motif from 3' end, bp") +
    ylab("Count") +
    labs(title="Distance of Detected PolyA Motif From 3' end") +
    mytheme+
    theme(legend.title=element_blank())
p.polyA_dist_subcat <- ggplot(data.FSMISM, aes(x=polyA_dist, color=subcategory)) +
  geom_freqpoly(binwidth=1, size=1) +
  scale_color_manual(values = subcat.palette)+
  xlab("Distance of polyA motif from 3' end, bp") +
  ylab("Count") +
  labs(title="Distance of Detected PolyA Motif From 3'End \n\nby FSM and ISM Subcategories")+
  mytheme+
  theme(legend.title=element_blank())

data.s2=rbind(data.other, data.NICNNC)
p.polyA_dist_subcat.s2 <- ggplot(data.s2, aes(x=polyA_dist, color=subcategory)) +
  geom_freqpoly(binwidth=1, size=1) +
  scale_color_manual(values = subcat.palette)+
  xlab("Distance of polyA motif from 3' end, bp") +
  ylab("Count") +
  labs(title="Distance of Detected PolyA Motif From 3'End \n\nby Non-FSM/ISM  Subcategories")+
  mytheme+
  theme(legend.title=element_blank())
}

### Bad quality control attributes
if (nrow(data.junction) > 0){
  t3.data.sets <- list()
  t3.list <- list()
  # (Fran) ToDo: USE COVERAGE DATA LATER
  # for FSM, ISM, NIC, and NNC, plot the percentage of RTS and non-canonical junction
  x <- filter(data.class, structural_category %in% c("FSM", "ISM", "NIC", "NNC" ) & exons > 1)
  
  t1.RTS <- group_by(x, structural_category, RTS_stage) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
  t2.RTS <- group_by(x, structural_category) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
  t3.RTS <- merge(t1.RTS, t2.RTS, by="structural_category")
  t3.RTS <- t3.RTS[-which(t3.RTS$structural_category=="ISM"),]
  t3.RTS$perc <- t3.RTS$count.x / t3.RTS$count.y * 100
  t3.RTS <- subset(t3.RTS, RTS_stage=='TRUE');
  n_t3.RTS <- dim(t3.RTS)[1];
  if (n_t3.RTS > 0) {
	  t3.RTS$Var <- "RT switching"
  }

  # Liz: this is a placeholder for dealing with all_canonical being NA instead of "Non-canonical"
  x[is.na(x$all_canonical), "all_canonical"] <- "Non-canonical"
  t1.SJ <- group_by(x, structural_category, all_canonical) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
  t3.SJ <- merge(t1.SJ, t2.RTS, by="structural_category")
  t3.SJ$perc <- t3.SJ$count.x / t3.SJ$count.y * 100
  t3.a.SJ <- subset(t3.SJ, all_canonical=='Canonical');
  t3.SJ <- subset(t3.SJ, all_canonical=='Non-canonical');
  n_t3.SJ <- dim(t3.SJ)[1];
  if (n_t3.SJ > 0) {
    t3.SJ$Var <- "Non-canonical"
    t3.a.SJ$Var <- 'Canonical'
  }
  

  if (!all(is.na(x$predicted_NMD))){
    x[which(x$predicted_NMD=="TRUE"),"predicted_NMD"]="Predicted NMD"
    x[which(x$predicted_NMD=="FALSE"),"predicted_NMD"]="Not NMD predicted"
    t1.NMD <- group_by(x, structural_category, predicted_NMD) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
    t3.NMD <- merge(t1.NMD, t2.RTS, by="structural_category")
    t3.NMD$perc <- t3.NMD$count.x / t3.NMD$count.y * 100
    t3.NMD <- subset(t3.NMD, predicted_NMD=='Predicted NMD');
    t3.NMD$Var=t3.NMD$predicted_NMD
  }
  if (!all(is.na(x$min_cov))){
    x[which(x$min_cov==0),"Coverage_SJ"]="Not Coverage SJ"
    x[which(x$min_cov>0),"Coverage_SJ"]="Coverage SJ"
    t1.Cov <- group_by(x, structural_category, Coverage_SJ) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
    t3.Cov <- merge(t1.Cov, t2.RTS, by="structural_category")
    t3.Cov$perc <- t3.Cov$count.x / t3.Cov$count.y * 100
    t3.a.Cov <- subset(t3.Cov, Coverage_SJ=='Coverage SJ');
    t3.Cov <- subset(t3.Cov, Coverage_SJ=='Not Coverage SJ');
    t3.Cov$Var=t3.Cov$Coverage_SJ
    t3.a.Cov$Var=t3.a.Cov$Coverage_SJ
    t3.data.sets[[length(t3.data.sets) + 1]]=x$min_cov
    t3.list[[length(t3.list) + 1]]=t3.a.Cov
  }
  
  if (!all(is.na(x$within_CAGE_peak))){
    x[which(!x$within_CAGE_peak),"Coverage_Cage"] <- "No Coverage CAGE"
    x[which(x$within_CAGE_peak),"Coverage_Cage"] <- "Has Coverage CAGE"
    t1.Cage <- group_by(x, structural_category, Coverage_Cage) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
    t3.Cage <- merge(t1.Cage, t2.RTS, by="structural_category")
    t3.Cage$perc <- t3.Cage$count.x / t3.Cage$count.y * 100
    t3.Cage <- subset(t3.Cage, Coverage_Cage=='Has Coverage CAGE');
    t3.Cage$Var <- t3.Cage$Coverage_Cage
    t3.data.sets[[length(t3.data.sets) + 1]] <- data.class$dist_to_cage_peak
    t3.list[[length(t3.list) + 1]] <- t3.Cage
    p28.a.Cage <- ggplot(t3.Cage, aes(x=structural_category, y=perc)) +
      geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[2] ,color="black") +
      geom_text(label=paste(round(t3.Cage$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) + 
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Isoforms, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle("CAGE Support\n\n") +
      theme(legend.title = element_blank())
  }
  if (!all(is.na(data.class$polyA_motif))) {
    x[which(is.na(x$polyA_motif)),"Coverage_PolyA"] <- "No Coverage PolyA"
    x[which(!is.na(x$polyA_motif)),"Coverage_PolyA"] <- "Has Coverage PolyA"
    t1.PolyA <- group_by(x, structural_category, Coverage_PolyA) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
    t3.PolyA <- merge(t1.PolyA, t2.RTS, by="structural_category")
    t3.PolyA$perc <- t3.PolyA$count.x / t3.PolyA$count.y * 100
    t3.PolyA <- subset(t3.PolyA, Coverage_PolyA=='Has Coverage PolyA');
    t3.PolyA$Var <- t3.PolyA$Coverage_PolyA
    t3.data.sets[[length(t3.data.sets) + 1]]=data.class$polyA_motif
    t3.list[[length(t3.list) + 1]]=t3.PolyA
    p28.a.polyA <- ggplot(t3.PolyA, aes(x=structural_category, y=perc)) +
      geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[3] ,color="black") +
      geom_text(label=paste(round(t3.PolyA$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) + 
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Isoforms, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle("PolyA Support\n\n") +
      guides(fill = guide_legend(title = "QC Attributes") )
  }
  
  x[which(x$diff_to_gene_TSS<=50),"Annotation"] <- "Annotated"
  x[which(x$diff_to_gene_TSS>50),"Annotation"] <- "Not annotated"
  t1.annot <- group_by(x, structural_category, Annotation) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
  t3.annot <- merge(t1.annot, t2.RTS, by="structural_category")
  t3.annot$perc <- t3.annot$count.x / t3.annot$count.y * 100
  t3.annot <- subset(t3.annot, Annotation=='Annotated');
  t3.annot$Var=t3.annot$Annotation
  p28.a.annot <- ggplot(t3.annot, aes(x=structural_category, y=perc)) +
    geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[6] ,color="black") +
    geom_text(label=paste(round(t3.annot$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) + 
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    ylab("Isoforms, %") +
    xlab("") +
    mytheme +
    theme(legend.position="bottom", axis.title.x = element_blank()) +
    ggtitle("Annotation Support\n\n") 
  #p28.Cov <- ggplot(t3.Cov, aes(x=structural_category, y=perc)) +
  #  geom_col(position='dodge', width = 0.7,  size=0.3, fill='lightblue', color="black") +
  #  geom_text(label=paste(round(t3.SJ$perc, 1),"%",sep=''), nudge_y=0.5) +
  #  scale_fill_manual(values = myPalette[11]) +
  #  ylab("Isoforms, %") +
  #  xlab("") +
  #  mytheme +
  #  theme(legend.position="bottom", axis.title.x = element_blank()) +
  #  ggtitle("Incidence of Non-Canonical Junctions\n\n") +
  #  guides(fill = guide_legend(title = "QC Attributes") )

  p28.RTS <- ggplot(t3.RTS, aes(x=structural_category, y=perc)) +
    geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[11], color="black") +
    geom_text(label=paste(round(t3.RTS$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
    scale_fill_manual(values = myPalette[9:11]) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    ylab("Isoforms, %") +
    xlab("") +
    mytheme +
    theme(legend.position="bottom", axis.title.x = element_blank()) +
    ggtitle("RT-switching\n\n")

  p28.SJ <- ggplot(t3.SJ, aes(x=structural_category, y=perc)) +
    geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[9] ,color="black") +
    geom_text(label=paste(round(t3.SJ$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) + 
    scale_fill_manual(values = myPalette[9:11]) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    ylab("Isoforms, %") +
    xlab("") +
    mytheme +
    theme(legend.position="bottom", axis.title.x = element_blank()) +
    ggtitle("Non-Canonical Junctions\n\n")
  p28.a.SJ <- ggplot(t3.a.SJ, aes(x=structural_category, y=perc)) +
    geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[7] ,color="black") +
    geom_text(label=paste(round(t3.a.SJ$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) + 
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    ylab("Isoforms, %") +
    xlab("") +
    mytheme +
    theme(legend.position="bottom", axis.title.x = element_blank()) +
    ggtitle("All Canonical Junctions\n\n")

  if (n_t3.SJ>0 & n_t3.RTS>0 & !all(is.na(x$min_cov)) & all(is.na(x$predicted_NMD))){
    p28.Cov <- ggplot(t3.Cov, aes(x=structural_category, y=perc)) +
      geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[10], color="black") +
      geom_text(label=paste(round(t3.Cov$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
      scale_fill_manual(values = myPalette[9:11]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Isoforms, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle("Splice Junctions Without Short Reads Coverage\n\n")
    p28.a.Cov <- ggplot(t3.a.Cov, aes(x=structural_category, y=perc)) +
      geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[10], color="black") +
      geom_text(label=paste(round(t3.a.Cov$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
      scale_fill_manual(values = myPalette[9:11]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Isoforms, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle("Splice Junctions With Short Reads Coverage\n\n") 


    t3 <- rbind(t3.RTS[,c(1,5,6)],t3.SJ[,c(1,5,6)], t3.Cov[,c(1,5,6)])

    p28 <- ggplot(data=t3, aes(x=structural_category, y=perc, fill= Var)) +
      geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
      scale_fill_manual(values = myPalette[9:11]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Transcripts, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle( "Summary Features of Bad Quality\n\n" ) +
      theme(legend.title = element_blank())
    
    #good quality control
    t3.a <- rbind(t3.annot[,c(1,5,6)], t3.a.SJ[,c(1,5,6)], t3.a.Cov[,c(1,5,6)])

  }else if (n_t3.SJ>0 & n_t3.RTS>0 & all(is.na(x$min_cov)) & all(is.na(x$predicted_NMD))) {
    t3=rbind(t3.RTS[,c(1,5,6)],t3.SJ[,c(1,5,6)])
    p28 <- ggplot(data=t3, aes(x=structural_category, y=perc, fill= Var)) +
      geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
      scale_fill_manual(values = myPalette[c(9,11)]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Transcripts, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle( "Quality Control Attributes Across Structural Categories\n\n" ) +
      theme(legend.title = element_blank())
    #good quality control
    t3.a=rbind(t3.annot[,c(1,5,6)], t3.a.SJ[,c(1,5,6)])

  }else if (n_t3.SJ>0 & n_t3.RTS>0 & all(is.na(x$min_cov)) & !all(is.na(x$predicted_NMD))){
    p28.NMD <- ggplot(t3.NMD, aes(x=structural_category, y=perc)) +
      geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[5], color="black") +
      geom_text(label=paste(round(t3.NMD$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
      scale_fill_manual(values = myPalette[9:11]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Isoforms, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle("Nonsense-Mediated Decay by Structural Category\n\n")
    t3=rbind(t3.RTS[,c(1,5,6)],t3.SJ[,c(1,5,6)], t3.NMD[,c(1,5,6)])
    p28 <- ggplot(data=t3, aes(x=structural_category, y=perc, fill= Var)) +
      geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
      scale_fill_manual(values = myPalette[c(9,5,11)]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Transcripts, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle( "Quality Control Attributes Across Structural Categories\n\n" ) +
      theme(legend.title = element_blank())
    #good quality control
    t3.a=rbind(t3.annot[,c(1,5,6)], t3.a.SJ[,c(1,5,6)])
  
  }else if (n_t3.SJ>0 & n_t3.RTS>0) {
    p28.NMD <- ggplot(t3.NMD, aes(x=structural_category, y=perc)) +
      geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[5], color="black") +
      geom_text(label=paste(round(t3.NMD$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
      scale_fill_manual(values = myPalette[9:11]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Isoforms, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle("Nonsense-Mediated Decay by Structural Category\n\n")
    p28.Cov <- ggplot(t3.Cov, aes(x=structural_category, y=perc)) +
      geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[10], color="black") +
      geom_text(label=paste(round(t3.Cov$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
      scale_fill_manual(values = myPalette[9:11]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Isoforms, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle("Splice Junctions Without Short Read Coverage\n\n")
    p28.a.Cov <- ggplot(t3.a.Cov, aes(x=structural_category, y=perc)) +
      geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[10], color="black") +
      geom_text(label=paste(round(t3.a.Cov$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
      scale_fill_manual(values = myPalette[9:11]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Isoforms, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle("Splice Junctions With Short Read Coverage\n\n")
    t3=rbind(t3.RTS[,c(1,5,6)],t3.SJ[,c(1,5,6)],t3.Cov[,c(1,5,6)], t3.NMD[,c(1,5,6)])
    p28 <- ggplot(data=t3, aes(x=structural_category, y=perc, fill= Var)) +
      geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
      scale_fill_manual(values = myPalette[c(9,10,5,11)]) +
      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
      ylab("Transcripts, %") +
      xlab("") +
      mytheme +
      theme(legend.position="bottom", axis.title.x = element_blank()) +
      ggtitle( "Quality Control Attributes Across Structural Categories\n\n" ) +
      theme(legend.title = element_blank())
    #good quality control
    t3.a <- rbind(t3.annot[,c(1,5,6)], t3.a.SJ[,c(1,5,6)], t3.a.Cov[,c(1,5,6)])
    
  }

}
t3.aa <-  rbind(t3.annot[,c("structural_category", "perc", "Var")], t3.a.SJ[,c(1,5,6)])

for(i in 1:length(t3.list)){
  set=data.frame(t3.data.sets[i])
  c=data.frame(t3.list[i])
  if (!all(is.na(set))){
    t.temp=t3.aa
    t3.aa = rbind(t.temp, c[,c(1,5,6)])
  }
}


p28.a <- ggplot(data=t3.aa, aes(x=structural_category, y=perc, fill= Var)) +
  geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  scale_fill_manual(values = c(myPalette)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
  ylab("Transcripts, %") +
  xlab("") +
  mytheme +
  theme(legend.position="bottom", axis.title.x = element_blank()) +
  ggtitle( "Good Quality Control Attributes Across Structural Categories\n\n" ) +
  theme(axis.text.y = element_text(size=10),
        axis.text.x  = element_text(size=10))+
  theme(legend.title = element_blank())
 


#TSS ratio
data.ratio = rbind(data.refmatch[,c(6,47)], data.ISM[,c(6,47)])
if (!all(is.na(data.ratio$ratio_TSS))){
  require(scales)
  p28.a.ratio=ggplot(data.ratio, aes(x=ratio_TSS, fill=structural_category)) + 
    geom_density(alpha=0.6)+
    labs(x="TSS ratio, log2", y="Density", title="TSS Ratio\n\nFSM Reference Match vs ISM\n\n") +
    scale_fill_manual(values = myPalette, breaks=c("FSM", "ISM"),
                      labels=c("FSM reference match", "ISM"), drop=F)+
    geom_vline(xintercept=1, linetype="dashed", color = "red")+
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    mytheme +
    theme(legend.position="bottom", legend.title = element_blank())
  p28.a.ratio.pdf <- p28.a.ratio + 
    scale_x_continuous(trans='log2', breaks = trans_breaks("log2", function(x) 2^x),
                                                      labels = trans_format("log2", math_format(2^.x)))
  p28.a.ratio.html <- p28.a.ratio + 
    scale_x_continuous(trans='log2')
}

# PLOT p30,p31,p32: percA by subcategory


p30.s1 <- ggplot(data=data.FSMISM, 
                 aes(y=perc_A_downstream_TTS, x=subcategory, fill=subcategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.1, position = "dodge") +
  mytheme +
  xlab("Structural category") +
  ylab("'A's, %") +
  labs(title="Possible Intra-Priming by Structural Category\n\n",
       subtitle="Percent of genomic 'A's in downstream 20 bp\n\n") +
  theme(legend.position="right", legend.title=element_blank()) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(strip.background = element_rect(color = "white"),
        strip.placement = "outside", strip.text = element_text(size = 10)) +
  scale_fill_manual(values=subcat.palette, drop=T) +
  facet_grid(~structural_category, scales = "free", switch = "x")

p30.s2 <- ggplot(data=data.NICNNC, 
                 aes(y=perc_A_downstream_TTS, x=subcategory, fill=subcategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.1, position = "dodge") +
  mytheme +
  xlab("Structural category") +
  ylab("'A's, %") +
  labs(title="Possible Intra-Priming by Structural Category\n\n",
       subtitle="Percent of genomic 'A's in downstream 20 bp\n\n") +
  theme(legend.position="right", legend.title=element_blank()) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(strip.background = element_rect(color = "white"),
        strip.placement = "outside", strip.text = element_text(size = 10)) +
  scale_fill_manual(values=subcat.palette, drop=T) +
  facet_grid(~structural_category, scales = "free", switch = "x")

p30.s3 <- ggplot(data=data.other, 
                 aes(y=perc_A_downstream_TTS, x=subcategory, fill=subcategory)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.1, position = "dodge") +
  mytheme +
  xlab("Structural category") +
  ylab("'A's, %") +
  labs(title="Possible Intra-Priming by Structural Category\n\n",
       subtitle="Percent of genomic 'A's in downstream 20 bp\n\n") +
  theme(legend.position="right", legend.title=element_blank()) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(strip.background = element_rect(color = "white"),
        strip.placement = "outside", strip.text = element_text(size = 10)) +
  scale_fill_manual(values=subcat.palette, drop=T) +
  facet_grid(~structural_category, scales = "free", switch = "x")
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

p31 <- ggplot(data=data.class, aes(y=perc_A_downstream_TTS, x=structural_category, fill=exonCat)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
  mytheme +
  scale_fill_manual(breaks=c("Mono-Exon", "Multi-Exon"),
                    labels=c("Mono-Exon Isoforms", "Multi-Exon Isoforms"), values=myPalette) +
  ylab("'A's, %") +
  theme(legend.position="bottom", legend.title=element_blank()) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
  labs(title = "Mono- vs Multi-Exon Possible Intra-Priming\n\n",
       subtitle = "Percent of genomic 'A's in downstream 20 bp\n\n") +
  theme(axis.title.x=element_blank())

p32 <- ggplot(data=data.class, aes(y=perc_A_downstream_TTS, x=structural_category, fill=coding)) +
  geom_boxplot(color="black", size=0.3, outlier.size = 0.2) + mytheme +
  scale_fill_manual(breaks=c("Coding", "Non coding"),
                    labels=c("Coding Isoforms", "Non-Coding Isoforms"), values=myPalette[3:4]) +
  ylab("'A's, % ") +
  theme(legend.position="bottom", legend.title=element_blank() ) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
  labs(title = "Coding vs Non-Coding Possible Intra-Priming\n\n",
       subtitle = "Percent of genomic 'A's in downstream 20 bp\n\n") +
  theme(axis.title.x=element_blank())

##### % of FSM or ISM associated to the same transcript ( histogram )

if (nrow(data.ISM) > 0 || nrow(data.FSM) > 0) {
  ism_per_transcript=data.ISM %>% group_by(associated_transcript, structural_category) %>% dplyr::summarize(dplyr::n())
  names(ism_per_transcript)[3]<-"ISM_per_tr"
  fsm_per_transcript=data.FSM %>% group_by(associated_transcript, structural_category) %>% dplyr::summarize(dplyr::n())
  names(fsm_per_transcript)[3]<-"FSM_per_tr"

  iso_per_knownTr=merge(x = fsm_per_transcript , y=ism_per_transcript, by = "associated_transcript", all=T)
  iso_per_knownTr$ISM_per_tr[is.na(iso_per_knownTr$ISM_per_tr)] <- 0
  iso_per_knownTr$FSM_per_tr[is.na(iso_per_knownTr$FSM_per_tr)] <- 0

  iso_per_knownTr$total_iso=apply(iso_per_knownTr, 1 , function(X) (as.integer(X[3]) + as.integer(X[5])) )
  #iso_per_knownTr$FSM_per_tr[is.na(iso_per_knownTr$FSM_per_tr)] <- 0
  #iso_per_knownTr$ISM_per_tr[is.na(iso_per_knownTr$ISM_per_tr)] <- 0
  #iso_per_knownTr$perc_FSM=apply(iso_per_knownTr, 1 , function(X) (as.numeric(X[3]) / as.numeric(X[6]) ) )
  #iso_per_knownTr$perc_FSM[is.na(iso_per_knownTr$perc_FSM)] <- 0
  #iso_per_knownTr$perc_FSM_range=cut(iso_per_knownTr$perc_FSM, breaks = c(-0.1,0.2,0.4,0.6,0.8,1),
  #                                   labels = c("FSM<=20% ", "20%<FSM<=40%", "40%<FSM<=60%", "60%<FSM<=80%","FSM>80%"))

  #
  iso_per_knownTr$FSM_cat=NA
  iso_per_knownTr$FSM_cat=apply(iso_per_knownTr,1, function(X){
     if(as.numeric(X["FSM_per_tr"])==1){
      return("Unique")
    }else if(as.numeric(X["FSM_per_tr"])>1){
      return("Multiple")
    }else{
      return("NULL")
    }})

  iso_per_knownTr$FSM_bin=apply(iso_per_knownTr,1, function(X){
    if(as.numeric(X["FSM_per_tr"])<8){
      return(as.character(X["FSM_per_tr"]))
    }else{
      return("8+")
    }})

  iso_per_knownTr$ISM_cat=NA
  iso_per_knownTr$ISM_cat=apply(iso_per_knownTr,1, function(X){
    if(as.numeric(X["ISM_per_tr"])==1){
      return("Unique")
    }else if(as.numeric(X["ISM_per_tr"])>1){
      return("Multiple")
    }else{
      return("NULL")
    }})

  iso_per_knownTr$ISM_bin=apply(iso_per_knownTr,1, function(X){
    if(as.numeric(X["ISM_per_tr"])<8){
      return(as.character(X["ISM_per_tr"]))
    }else{
      return("8+")
    }})


  iso_per_knownTr$total_cat=NA
  iso_per_knownTr$total_cat=apply(iso_per_knownTr,1, function(X){
    if(as.numeric(X["total_iso"])==1){
      return("Unique")
    }else if(as.numeric(X["total_iso"])>1){
      return("Multiple")
      }else{
        return("NULL")
      }})

  iso_per_knownTr$total_bin=apply(iso_per_knownTr,1, function(X){
    if(as.numeric(X["total_iso"])<8){
      return(as.character(X["total_iso"]))
    }else{
      return("8+")
    }})
  iso_per_knownTr$FSM_cat=factor(iso_per_knownTr$FSM_cat, levels=c("Unique", "Multiple"))
  iso_per_knownTr$ISM_cat=factor(iso_per_knownTr$ISM_cat, levels=c("Unique", "Multiple"))
  iso_per_knownTr$total_cat=factor(iso_per_knownTr$total_cat, levels=c("Unique", "Multiple"))

  max_y=max(table(iso_per_knownTr[which(iso_per_knownTr$FSM_cat!="NULL"),"FSM_cat"]))+10
  new.1.FSM <- ggplot(iso_per_knownTr[which(iso_per_knownTr$FSM_cat!="NULL"),])+
    geom_bar(aes(x=FSM_cat, fill=FSM_bin), color = "black", width = 0.5) +
    mytheme+
    geom_text(aes(x = FSM_cat, label = stat(count)), stat = "count", vjust = -0.5) +
    scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
    scale_fill_brewer("Total FSM \nper reference ID", palette = "Blues") +
    labs(x="FSM per reference transcript",
         y="Count of reference transcripts",
         title="Reference Transcript Redundancy",
         subtitle="Only FSM")

  max_y=max(table(iso_per_knownTr[which(iso_per_knownTr$ISM_cat!="NULL"),"ISM_cat"]))+10
  new.1.ISM <- ggplot(iso_per_knownTr[which(iso_per_knownTr$ISM_cat!="NULL"),])+
    geom_bar(aes(x=ISM_cat, fill=ISM_bin), color = "black", width = 0.5) +
    mytheme+
    geom_text(aes(x = ISM_cat, label = stat(count)), stat = "count", vjust = -0.5) +
    scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
    scale_fill_brewer("Total ISM \nper reference ID", palette = "Oranges") +
    labs(x="ISM per reference transcript",
         y="Count of Reference Transcripts",
         title="Reference Transcript Redundancy",
         subtitle="Only ISM")

  max_y=max(table(iso_per_knownTr[which(iso_per_knownTr$total_cat!="NULL"),"total_cat"]))+10
  new.1.total <- ggplot(iso_per_knownTr[which(iso_per_knownTr$total_cat!="NULL"),])+
    geom_bar(aes(x=total_cat, fill=total_bin), color = "black", width = 0.5) +
    mytheme+
    geom_text(aes(x = total_cat, label = stat(count)), stat = "count", vjust = -0.5) +
    scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
    scale_fill_brewer("Total FSM+ISM \nper reference", palette = "Greens") +
    labs(x="FSM+ISM per reference transcript",
         y="Count of reference transcripts",
         title="Reference Transcript Redundancy",
         subtitle="FSM+ISM")

  #new.1.p<-ggplot(iso_per_knownTr, aes(x=total_iso, fill=factor(perc_FSM_range))) +
  #  geom_histogram(binwidth = 1 , col="white") +
  #  scale_x_continuous(breaks = seq(0,max_num_iso,1)) +
  #  scale_fill_manual(values = myPalette, guide='none', name="Percentage of FSM isoforms") +
  #  mytheme+theme(axis.title.x=element_text()) +
  #  theme(legend.justification=c(1,1), legend.position=c(1,1))  +
  #  guides(fill = guide_legend(keywidth = 0.9, keyheight = 0.9)) +
  #  labs(x="Number of FSM+ISM isoforms \n associated to the same ref. transcript",
  #       y="Counts", title="Accumulation of FSM and ISM Isoforms \n\n Associated to the Same Reference Transcript.")

  #### Now only with cage + isoforms
  if (!all(is.na(data.class$dist_to_cage_peak))) {
    ism_per_transcript_cage=data.ISM[which(data.ISM$within_CAGE_peak),] %>% group_by(associated_transcript, structural_category) %>% dplyr::summarize(dplyr::n())
    names(ism_per_transcript_cage)[3]<-"ISM_per_tr"
    fsm_per_transcript_cage=data.FSM[which(data.FSM$within_CAGE_peak),] %>% group_by(associated_transcript, structural_category) %>% dplyr::summarize(dplyr::n())
    names(fsm_per_transcript_cage)[3]<-"FSM_per_tr"

    iso_per_knownTr_cage=merge(x = fsm_per_transcript_cage , y=ism_per_transcript_cage, by = "associated_transcript", all=T)
    iso_per_knownTr_cage$ISM_per_tr[is.na(iso_per_knownTr_cage$ISM_per_tr)] <- 0
    iso_per_knownTr_cage$FSM_per_tr[is.na(iso_per_knownTr_cage$FSM_per_tr)] <- 0

    iso_per_knownTr_cage$total_iso=apply(iso_per_knownTr_cage, 1 , function(X) as.integer(X[3]) + as.integer(X[5]) )
    #iso_per_knownTr_cage$perc_FSM=apply(iso_per_knownTr_cage, 1 , function(X) (as.numeric(X[3]) / as.numeric(X[6]) ) )
    #iso_per_knownTr_cage$perc_FSM[is.na(iso_per_knownTr_cage$perc_FSM)] <- 0
    #iso_per_knownTr_cage$perc_FSM_range=cut(iso_per_knownTr_cage$perc_FSM, breaks = c(-0.1,0.2,0.4,0.6,0.8,1),
    #                                        labels = c("FSM<=20% ", "20%<FSM<=40%", "40%<FSM<=60%", "60%<FSM<=80%","FSM>80%"))
    #max_num_iso_cage=max(iso_per_knownTr_cage$total_iso)

    iso_per_knownTr_cage$FSM_cat=NA
    iso_per_knownTr_cage$FSM_cat=apply(iso_per_knownTr_cage,1, function(X){
      if(as.numeric(X["FSM_per_tr"])==1){
        return("Unique")
      }else if(as.numeric(X["FSM_per_tr"])>1){
        return("Multiple")
      }else{
        return("NULL")
      }})

    iso_per_knownTr_cage$FSM_bin=apply(iso_per_knownTr_cage,1, function(X){
      if(as.numeric(X["FSM_per_tr"])<8){
        return(as.character(X["FSM_per_tr"]))
      }else{
        return("8+")
      }})

    iso_per_knownTr_cage$ISM_cat=NA
    iso_per_knownTr_cage$ISM_cat=apply(iso_per_knownTr_cage,1, function(X){
      if(as.numeric(X["ISM_per_tr"])==1){
        return("Unique")
      }else if(as.numeric(X["ISM_per_tr"])>1){
        return("Multiple")
      }else{
        return("NULL")
      }})

    iso_per_knownTr_cage$ISM_bin=apply(iso_per_knownTr_cage,1, function(X){
      if(as.numeric(X["ISM_per_tr"])<8){
        return(as.character(X["ISM_per_tr"]))
      }else{
        return("8+")
      }})


    iso_per_knownTr_cage$total_cat=NA
    iso_per_knownTr_cage$total_cat=apply(iso_per_knownTr_cage,1, function(X){
      if(as.numeric(X["total_iso"])==1){
        return("Unique")
      }else if(as.numeric(X["total_iso"])>1){
        return("Multiple")
      }else{
        return("NULL")
      }})

    iso_per_knownTr_cage$total_bin=apply(iso_per_knownTr_cage,1, function(X){
      if(as.numeric(X["total_iso"])<8){
        return(as.character(X["total_iso"]))
      }else{
        return("8+")
      }})

    iso_per_knownTr_cage$FSM_cat=factor(iso_per_knownTr_cage$FSM_cat, levels=c("Unique", "Multiple"))
    iso_per_knownTr_cage$ISM_cat=factor(iso_per_knownTr_cage$ISM_cat, levels=c("Unique", "Multiple"))
    iso_per_knownTr_cage$total_cat=factor(iso_per_knownTr_cage$total_cat, levels=c("Unique", "Multiple"))


    max_y=max(table(iso_per_knownTr_cage[which(iso_per_knownTr_cage$FSM_cat!="NULL"),"FSM_cat"]))+10
    new.2.FSM <- ggplot(iso_per_knownTr_cage[which(iso_per_knownTr_cage$FSM_cat!="NULL"),])+
      geom_bar(aes(x=FSM_cat, fill=FSM_bin), color = "black", width = 0.5) +
      mytheme+
      geom_text(aes(x = FSM_cat, label = stat(count)), stat = "count", vjust = -0.5) +
      scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
      scale_fill_brewer("Total FSM \nper reference ID", palette = "Blues") +
      labs(x="FSM per reference transcript",
           y="Count of reference transcripts",
           title="Reference Transcript Redundancy",
           subtitle="Only FSM with CAGE support")

    max_y=max(table(iso_per_knownTr_cage[which(iso_per_knownTr_cage$ISM_cat!="NULL"),"ISM_cat"]))+10
    new.2.ISM <- ggplot(iso_per_knownTr_cage[which(iso_per_knownTr_cage$ISM_cat!="NULL"),])+
      geom_bar(aes(x=ISM_cat, fill=ISM_bin), color = "black", width = 0.5) +
      mytheme+
      geom_text(aes(x = ISM_cat, label = stat(count)), stat = "count", vjust = -0.5) +
      scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
      scale_fill_brewer("Total ISM \nper reference ID", palette = "Oranges") +
      labs(x="ISM per reference transcript",
           y="Count of reference transcripts",
           title="Reference Transcript Redundancy",
           subtitle="Only ISM with CAGE support")

    max_y=max(table(iso_per_knownTr_cage[which(iso_per_knownTr_cage$total_cat!="NULL"),"total_cat"]))+10
    new.2.total <- ggplot(iso_per_knownTr_cage[which(iso_per_knownTr_cage$total_cat!="NULL"),])+
      geom_bar(aes(x=total_cat, fill=total_bin), color = "black", width = 0.5) +
      mytheme+
      geom_text(aes(x = total_cat, label = stat(count)), stat = "count", vjust = -0.5) +
      scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
      scale_fill_brewer("Total FSM+ISM \nper reference", palette = "Greens") +
      labs(x="FSM+ISM per reference transcript",
           y="Count of reference transcripts",
           title="Reference Transcript Redundancy",
           subtitle="FSM+ISM with CAGE support")

    #new.2.p<-ggplot(iso_per_knownTr_cage, aes(x=total_iso, fill=factor(perc_FSM_range))) +
    #  geom_histogram(binwidth = 1, col="white") +
    #  scale_x_continuous(breaks = seq(0,max_num_iso_cage,1)) +
    #  scale_fill_manual(values = myPalette, guide='none', name="Percentage of FSM isoforms") +
    #  mytheme+theme(axis.title.x=element_text()) +
    #  theme(legend.justification=c(1,1), legend.position=c(1,1))  +
    #  guides(fill = guide_legend(keywidth = 0.9, keyheight = 0.9)) +
    #  labs(x="Number of FSM+ISM isoforms \n associated to the same ref. transcript",
    #       y="Counts", title="Accumulation of FSM and ISM Isoforms \n\n Associated to the Same Reference Transcript.",
    #       subtitle="ONLY polyA motif and CAGE + isoforms")
  }

  ### Now only with polyA motif = T isoforms
  if (!all(is.na(data.class$polyA_motif))) {
    ism_per_transcript_polya=data.ISM[which(!is.na(data.ISM$polyA_motif)),] %>% group_by(associated_transcript, structural_category) %>% dplyr::summarize(dplyr::n(), .groups='keep')
    names(ism_per_transcript_polya)[3]<-"ISM_per_tr"
    fsm_per_transcript_polya=data.FSM[which(!is.na(data.FSM$polyA_motif)),] %>% group_by(associated_transcript, structural_category) %>% dplyr::summarize(dplyr::n(), .groups='keep')
    names(fsm_per_transcript_polya)[3]<-"FSM_per_tr"

    iso_per_knownTr_polya=merge(x = fsm_per_transcript_polya , y=ism_per_transcript_polya, by = "associated_transcript", all=T)
    iso_per_knownTr_polya$ISM_per_tr[is.na(iso_per_knownTr_polya$ISM_per_tr)] <- 0
    iso_per_knownTr_polya$FSM_per_tr[is.na(iso_per_knownTr_polya$FSM_per_tr)] <- 0
    iso_per_knownTr_polya$total_iso=apply(iso_per_knownTr_polya, 1 , function(X) as.integer(X[3]) + as.integer(X[5]) )
    #iso_per_knownTr_polya$perc_FSM=apply(iso_per_knownTr_polya, 1 , function(X) (as.numeric(X[3]) / as.numeric(X[6]) ) )
    #iso_per_knownTr_polya$perc_FSM[is.na(iso_per_knownTr_polya$perc_FSM)] <- 0
    #iso_per_knownTr_polya$perc_FSM_range=cut(iso_per_knownTr_polya$perc_FSM, breaks = c(-0.1,0.2,0.4,0.6,0.8,1),
    #                                         labels = c("FSM<=20% ", "20%<FSM<=40%", "40%<FSM<=60%", "60%<FSM<=80%","FSM>80%"))
    #max_num_iso_polya=max(iso_per_knownTr_polya$total_iso)

    iso_per_knownTr_polya$FSM_cat=NA
    iso_per_knownTr_polya$FSM_cat=apply(iso_per_knownTr_polya,1, function(X){
      if(as.numeric(X["FSM_per_tr"])==1){
        return("Unique")
      }else if(as.numeric(X["FSM_per_tr"])>1){
        return("Multiple")
      }else{
        return("NULL")
      }})

    iso_per_knownTr_polya$FSM_bin=apply(iso_per_knownTr_polya,1, function(X){
      if(as.numeric(X["FSM_per_tr"])<8){
        return(as.character(X["FSM_per_tr"]))
      }else{
        return("8+")
      }})

    iso_per_knownTr_polya$ISM_cat=NA
    iso_per_knownTr_polya$ISM_cat=apply(iso_per_knownTr_polya,1, function(X){
      if(as.numeric(X["ISM_per_tr"])==1){
        return("Unique")
      }else if(as.numeric(X["ISM_per_tr"])>1){
        return("Multiple")
      }else{
        return("NULL")
      }})

    iso_per_knownTr_polya$ISM_bin=apply(iso_per_knownTr_polya,1, function(X){
      if(as.numeric(X["ISM_per_tr"])<8){
        return(as.character(X["ISM_per_tr"]))
      }else{
        return("8+")
      }})


    iso_per_knownTr_polya$total_cat=NA
    iso_per_knownTr_polya$total_cat=apply(iso_per_knownTr_polya,1, function(X){
      if(as.numeric(X["total_iso"])==1){
        return("Unique")
      }else if(as.numeric(X["total_iso"])>1){
        return("Multiple")
      }else{
        return("NULL")
      }})

    iso_per_knownTr_polya$total_bin=apply(iso_per_knownTr_polya,1, function(X){
      if(as.numeric(X["total_iso"])<8){
        return(as.character(X["total_iso"]))
      }else{
        return("8+")
      }})

    iso_per_knownTr_polya$FSM_cat=factor(iso_per_knownTr_polya$FSM_cat, levels=c("Unique", "Multiple"))
    iso_per_knownTr_polya$ISM_cat=factor(iso_per_knownTr_polya$ISM_cat, levels=c("Unique", "Multiple"))
    iso_per_knownTr_polya$total_cat=factor(iso_per_knownTr_polya$total_cat, levels=c("Unique", "Multiple"))

    max_y=max(table(iso_per_knownTr_polya[which(iso_per_knownTr_polya$FSM_cat!="NULL"),"FSM_cat"]))+10
    new.3.FSM <- ggplot(iso_per_knownTr_polya[which(iso_per_knownTr_polya$FSM_cat!="NULL"),])+
      geom_bar(aes(x=FSM_cat, fill=FSM_bin), color = "black", width = 0.5) +
      mytheme+
      geom_text(aes(x = FSM_cat, label = stat(count)), stat = "count", vjust = -0.5) +
      scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
      scale_fill_brewer("Total FSM \nper reference ID", palette = "Blues") +
      labs(x="FSM per reference transcript",
           y="Count of reference transcripts",
           title="Reference Transcript Redundancy",
           subtitle="Only FSM with a polyA motif found")

    max_y=max(table(iso_per_knownTr_polya[which(iso_per_knownTr_polya$ISM_cat!="NULL"),"ISM_cat"]))+10
    new.3.ISM <- ggplot(iso_per_knownTr_polya[which(iso_per_knownTr_polya$ISM_cat!="NULL"),])+
      geom_bar(aes(x=ISM_cat, fill=ISM_bin), color = "black", width = 0.5) +
      mytheme+
      geom_text(aes(x = ISM_cat, label = stat(count)), stat = "count", vjust = -0.5) +
      scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
      scale_fill_brewer("Total ISM \nper reference ID", palette = "Oranges") +
      labs(x="ISM per reference transcript",
           y="Count of reference transcripts",
           title="Reference Transcript Redundancy",
           subtitle="Only ISM with a polyA motif found")

    max_y=max(table(iso_per_knownTr_polya[which(iso_per_knownTr_polya$total_cat!="NULL"),"total_cat"]))+10
    new.3.total <- ggplot(iso_per_knownTr_polya[which(iso_per_knownTr_polya$total_cat!="NULL"),])+
      geom_bar(aes(x=total_cat, fill=total_bin), color = "black", width = 0.5) +
      mytheme+
      geom_text(aes(x = total_cat, label = stat(count)), stat = "count", vjust = -0.5) +
      scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
      scale_fill_brewer("Total FSM+ISM \nper reference", palette = "Greens") +
      labs(x="FSM+ISM per reference transcript",
           y="Count of reference transcripts",
           title="Reference Transcript Redundancy",
           subtitle="FSM+ISM with a polyA motif found")

    #new.3.p<-ggplot(iso_per_knownTr_polya, aes(x=total_iso, fill=factor(perc_FSM_range))) +
    #  geom_histogram(binwidth = 1, col="white") +
    #  scale_x_continuous(breaks = seq(0,max_num_iso_polya,1)) +
    #  scale_fill_manual(values = myPalette, guide='none', name="Percentage of FSM isoforms") +
    #  mytheme+theme(axis.title.x=element_text()) +
    #  theme(legend.justification=c(1,1), legend.position=c(1,1))  +
    #  guides(fill = guide_legend(keywidth = 0.9, keyheight = 0.9)) +
    #  labs(x="Number of FSM+ISM isoforms \n associated to the same ref. transcript",
    #       y="Counts", title="Accumulation of FSM and ISM Isoforms \n\n Associated to the Same Reference Transcript.",
    #       subtitle="ONLY polyA motif and CAGE + isoforms")

  }


  #### Now with just isoforms polyA and Cage +
  if (!all(is.na(data.class$polyA_motif)) && !all(is.na(data.class$dist_to_cage_peak))) {
    ism_per_transcript_cage_polya=data.ISM[which(!is.na(data.ISM$polyA_motif) & data.ISM$within_CAGE_peak),] %>% group_by(associated_transcript, structural_category) %>% dplyr::summarize(dplyr::n(), .groups='keep')
    names(ism_per_transcript_cage_polya)[3]<-"ISM_per_tr"
    fsm_per_transcript_cage_polya=data.FSM[which(!is.na(data.FSM$polyA_motif) & data.FSM$within_CAGE_peak),] %>% group_by(associated_transcript, structural_category) %>% dplyr::summarize(dplyr::n(), .groups='keep')
    names(fsm_per_transcript_cage_polya)[3]<-"FSM_per_tr"

    iso_per_knownTr_cage_polya=merge(x = fsm_per_transcript_cage_polya , y=ism_per_transcript_cage_polya, by = "associated_transcript", all=T)
    iso_per_knownTr_cage_polya$ISM_per_tr[is.na(iso_per_knownTr_cage_polya$ISM_per_tr)] <- 0
    iso_per_knownTr_cage_polya$FSM_per_tr[is.na(iso_per_knownTr_cage_polya$FSM_per_tr)] <- 0
    iso_per_knownTr_cage_polya$total_iso=apply(iso_per_knownTr_cage_polya, 1 , function(X) as.integer(X[3]) + as.integer(X[5]) )
    #iso_per_knownTr_cage_polya$perc_FSM=apply(iso_per_knownTr_cage_polya, 1 , function(X) (as.numeric(X[3]) / as.numeric(X[6]) ) )
    #iso_per_knownTr_cage_polya$perc_FSM[is.na(iso_per_knownTr_cage_polya$perc_FSM)] <- 0
    #iso_per_knownTr_cage_polya$perc_FSM_range=cut(iso_per_knownTr_cage_polya$perc_FSM, breaks = c(-0.1,0.2,0.4,0.6,0.8,1),
    #                                              labels = c("FSM<=20% ", "20%<FSM<=40%", "40%<FSM<=60%", "60%<FSM<=80%","FSM>80%"))
    #max_num_iso_cage_polya=max(iso_per_knownTr_cage_polya$total_iso)

    iso_per_knownTr_cage_polya$FSM_cat=NA
    iso_per_knownTr_cage_polya$FSM_cat=apply(iso_per_knownTr_cage_polya,1, function(X){
      if(as.numeric(X["FSM_per_tr"])==1){
        return("Unique")
      }else if(as.numeric(X["FSM_per_tr"])>1){
        return("Multiple")
      }else{
        return("NULL")
      }})

    iso_per_knownTr_cage_polya$FSM_bin=apply(iso_per_knownTr_cage_polya,1, function(X){
      if(as.numeric(X["FSM_per_tr"])<8){
        return(as.character(X["FSM_per_tr"]))
      }else{
        return("8+")
      }})

    iso_per_knownTr_cage_polya$ISM_cat=NA
    iso_per_knownTr_cage_polya$ISM_cat=apply(iso_per_knownTr_cage_polya,1, function(X){
      if(as.numeric(X["ISM_per_tr"])==1){
        return("Unique")
      }else if(as.numeric(X["ISM_per_tr"])>1){
        return("Multiple")
      }else{
        return("NULL")
      }})

    iso_per_knownTr_cage_polya$ISM_bin=apply(iso_per_knownTr_cage_polya,1, function(X){
      if(as.numeric(X["ISM_per_tr"])<8){
        return(as.character(X["ISM_per_tr"]))
      }else{
        return("8+")
      }})


    iso_per_knownTr_cage_polya$total_cat=NA
    iso_per_knownTr_cage_polya$total_cat=apply(iso_per_knownTr_cage_polya,1, function(X){
      if(as.numeric(X["total_iso"])==1){
        return("Unique")
      }else if(as.numeric(X["total_iso"])>1){
        return("Multiple")
      }else{
        return("NULL")
      }})

    iso_per_knownTr_cage_polya$total_bin=apply(iso_per_knownTr_cage_polya,1, function(X){
      if(as.numeric(X["total_iso"])<8){
        return(as.character(X["total_iso"]))
      }else{
        return("8+")
      }})

    iso_per_knownTr_cage_polya$FSM_cat=factor(iso_per_knownTr_cage_polya$FSM_cat, levels=c("Unique", "Multiple"))
    iso_per_knownTr_cage_polya$ISM_cat=factor(iso_per_knownTr_cage_polya$ISM_cat, levels=c("Unique", "Multiple"))
    iso_per_knownTr_cage_polya$total_cat=factor(iso_per_knownTr_cage_polya$total_cat, levels=c("Unique", "Multiple"))

    max_y=max(table(iso_per_knownTr_cage_polya[which(iso_per_knownTr_cage_polya$FSM_cat!="NULL"),"FSM_cat"]))+10
    new.4.FSM <- ggplot(iso_per_knownTr_cage_polya[which(iso_per_knownTr_cage_polya$FSM_cat!="NULL"),])+
      geom_bar(aes(x=FSM_cat, fill=FSM_bin), color = "black", width = 0.5) +
      mytheme+
      geom_text(aes(x = FSM_cat, label = stat(count)), stat = "count", vjust = -0.5) +
      scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
      scale_fill_brewer("Total FSM \nper reference ID", palette = "Blues") +
      labs(x="FSM per reference transcript",
           y="Count of reference transcripts",
           title="Reference Transcript Redundancy",
           subtitle="Only FSM with CAGE support and polyA motif")

    max_y=max(table(iso_per_knownTr_cage_polya[which(iso_per_knownTr_cage_polya$ISM_cat!="NULL"),"ISM_cat"]))+10
    new.4.ISM <- ggplot(iso_per_knownTr_cage_polya[which(iso_per_knownTr_cage_polya$ISM_cat!="NULL"),])+
      geom_bar(aes(x=ISM_cat, fill=ISM_bin), color = "black", width = 0.5) +
      mytheme+
      geom_text(aes(x = ISM_cat, label = stat(count)), stat = "count", vjust = -0.5) +
      scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
      scale_fill_brewer("Total ISM \nper reference ID", palette = "Oranges") +
      labs(x="ISM per reference transcript",
           y="Count of reference transcripts",
           title="Reference Transcript Redundancy",
           subtitle="Only ISM with CAGE support and polyA motif")

    max_y=max(table(iso_per_knownTr_cage_polya[which(iso_per_knownTr_cage_polya$total_cat!="NULL"),"total_cat"]))+10
    new.4.total <- ggplot(iso_per_knownTr_cage_polya[which(iso_per_knownTr_cage_polya$total_cat!="NULL"),])+
      geom_bar(aes(x=total_cat, fill=total_bin), color = "black", width = 0.5) +
      mytheme+
      geom_text(aes(x = total_cat, label = stat(count)), stat = "count", vjust = -0.5) +
      scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
      scale_fill_brewer("Total FSM+ISM \nper reference", palette = "Greens") +
      labs(x="FSM+ISM per reference transcript",
           y="Count of reference transcripts",
           title="Reference Transcript Redundancy",
           subtitle="FSM+ISM with CAGE support and polyA motif")

    #new.4.p<-ggplot(iso_per_knownTr_cage_polya, aes(x=total_iso, fill=factor(perc_FSM_range))) +
    #  geom_histogram(binwidth = 1, col="white") +
    #  scale_x_continuous(breaks = seq(0,max_num_iso_cage_polya,1)) +
    #  scale_fill_manual(values = myPalette, guide='none', name="Percentage of FSM isoforms") +
    #  mytheme+theme(axis.title.x=element_text()) +
    #  theme(legend.justification=c(1,1), legend.position=c(1,1))  +
    #  guides(fill = guide_legend(keywidth = 0.9, keyheight = 0.9)) +
    #  labs(x="Number of FSM+ISM isoforms \n associated to the same Ref. transcript",
    #       y="Counts", title="Accumulation of FSM and ISM Isoforms \n\n Associated to the Same Reference Transcript.",
    #       subtitle="ONLY polyA motif and CAGE + isoforms")

  }

}



#### Rarefraction plots
if (saturation.curves=='True'){
  if (!all(is.na(data.class$FL))) {
    FL.counts <- as.matrix(data.class$FL)
    rownames(FL.counts) <- data.class$isoform
    colnames(FL.counts) <- "FL"
    myfactors <- data.frame(sample = c("FL"))
    rownames(myfactors) = colnames(FL.counts)
    mybiotype=as.matrix(data.class$coding)
    rownames(mybiotype)=data.class$isoform
    mycategory=as.matrix(data.class$structural_category)
    rownames(mycategory)=data.class$isoform
    mydata = readData(data = FL.counts, factors = myfactors, biotype = mybiotype, category=mycategory)

    rarefact <- LR.rarefaction(mydata , samples = 1)
    rar1 <- suppressWarnings(plot.rarefaction(rarefact, sample = 1, k = 1, depth.increase = 2, break.by = "category"))
    rar2 <- suppressWarnings(plot.rarefaction(rarefact, sample = 1, k = 2, depth.increase = 2, break.by = "category"))
    rar3 <- suppressWarnings(plot.rarefaction(rarefact, sample = 1, k = 3, depth.increase = 2, break.by = "category"))
    rar5 <- suppressWarnings(plot.rarefaction(rarefact, sample = 1, k = 5, depth.increase = 2, break.by = "category"))
  }
}



# PLOT pn1.2: Splice Junction relative coverage (if coverage and expression provided)
##### NEEDS TRANSCRIPT_COORD VALUES IN JUNCTIONS FILE (???)

#if (nrow(data.junction) > 0){
#  if (!all(is.na(data.junction$total_coverage)) & !all(is.na(data.class$iso_exp))){

#   data.junction$isoExp = data.class[data.junction$isoform, "iso_exp"]

#    total = aggregate(cbind(total_coverage,isoExp,transcript_coord) ~ junctionLabel, data = data.junction,
#                      FUN = function(x) c(mn = sum(x), n = min(x) ) )

#    total$relCov = total$total_coverage[,"n"] / total$isoExp[,"mn"]
#    total$minTSS = total$transcript_coord[,"n"]

#    uniqJunc = unique(data.junction[,c("junctionLabel", "canonical_known", "total_coverage")])
#    uniqJunc$notCov = uniqJunc$total_coverage == 0

#    uniqueJunc_nonCov = as.data.frame(table(uniqJunc[uniqJunc$totalCoverage==0,"canonical_known"])/table(uniqJunc$canonical_known)*100)

#    uniqJunc2 = merge(total, uniqJunc, by=1)
#    uniqJunc2$TSSrange =cut(uniqJunc2$minTSS, breaks = c(0,40,80,120,160,200,10000000), labels = c("0-40", "41-80", "81-120", "121-160", "161-200",">200"))



    # calculate total expression associated to each unique junction
#    sumExpPerJunc = tapply(data.junction$isoExp, data.junction$junctionLabel, sum)

#    data.junction$sumIsoExp = sumExpPerJunc[data.junction$junctionLabel]

#    data.junction$relCov = data.junction$total_coverage / data.junction$sumIsoExp

#    max_dist = max(data.junction$transcript_coord) +1

#    data.junction$TSSrange = cut(data.junction$transcript_coord, breaks = c(0,20,40,60,80,100,120,140,160,180,200,max_dist), labels = c("0-20", "21-40","41-80","61-80", "81-100","101-120", "121-140","141-160", "161-180", "181-200", ">200"))

#    pn1.2 <-ggplot(data=data.junction[data.junction$relCov<1,], aes(y=relCov,x=TSSrange,fill=canonical_known)) +
#      geom_boxplot(outlier.size = 0.2, size=0.3) +
#      scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
#      ylab("Relative coverage") +
#      xlab("# TSS distance range") +
#      mytheme_bw +
#      theme(legend.position="bottom", legend.title=element_blank())  +
#      ggtitle( "Junctions Relative Coverage\n\n\n") +
#      theme(axis.text.x = element_text(angle = 45,margin=margin(15,0,0,0), size=12))


#  }else{    uniqJunc = unique(data.junction[,c("junctionLabel", "canonical_known")])
#  }
#}

# table data arrange 
if (sum(!is.na(data.class$polyA_dist)) > 10) {
  df.polyA <- as.data.frame(group_by(data.class, by=structural_category) %>%
                              dplyr::summarise(count=dplyr::n(),
                                               polyA_detected=sum(!is.na(polyA_motif)),
                                               polyA_detected_perc=round(polyA_detected*100/count) ,
                                               .groups = 'keep'))
  df.polyA_freq <- as.data.frame(sort(table(data.class$polyA_motif),decreasing=T))
  df.polyA_freq$perc <- round(df.polyA_freq$Freq*100/sum(df.polyA_freq$Freq),1)
  
  df.polyA_subcat <- as.data.frame(group_by(data.class, by=subcategory) %>%
                                     dplyr::summarise(count=dplyr::n(),
                                                      polyA_detected=sum(!is.na(polyA_motif)),
                                                      polyA_detected_perc=round(polyA_detected*100/count) ,
                                                      .groups = 'keep'))
}




if (!all(is.na(data.class$dist_to_CAGE_peak))){
  df.cage <- as.data.frame(group_by(data.class, by=structural_category) %>%
                             dplyr::summarise(count=dplyr::n(),
                                              cage_detected=length(which(within_CAGE_peak)),
                                              cage_detected_perc=round(cage_detected*100/count) , 
                                              .groups = 'keep'))
  
  df.cage_subc <- as.data.frame(group_by(data.class, by=subcategory) %>%
                                  dplyr::summarise(count=dplyr::n(),
                                                   cage_detected=length(which(within_CAGE_peak)),
                                                   cage_detected_perc=round(cage_detected*100/count) ,
                                                   .groups = 'keep'))
}

###** Output plots

if (report.format == 'both') {
  invisible(generatePDFreport())
  rmarkdown::render(input = paste(utilities.path, "/report_qc/SQANTI3_report.Rmd", sep = "/"), 
                    intermediates_dir = output_directory, 
                    output_dir =  output_directory, 
                    output_file=html.report.file)
} else if (args[6] == 'pdf' & args[6] != 'html'){
  invisible(generatePDFreport())
} else {
  rmarkdown::render(input = paste(utilities.path, "/report_qc/SQANTI3_report.Rmd", sep = "/"), 
                    intermediates_dir = output_directory,
                    output_dir =  output_directory, 
                    output_file=html.report.file)
}

