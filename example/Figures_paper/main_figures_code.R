library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(RColorConesa)
library(scales)
library(patchwork)
library(ggpubr)
library(gghalves)
library(stringr)
library(MetBrewer)

if (!file.exists("main_figures")) {
  dir.create("main_figures")
}

cat.palette = c( "FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                 "NNC"="#EE6A50", "GenicGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                 "Intergenic" = "darksalmon", "GenicIntron"="#41B6C4")

subcat.palette = c("Alternative 3'end"='#02314d',
                   "Alternative 3'5'end"='#0e5a87',
                   "Alterantive 5'end"='#7ccdfc',
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

pub_theme <- theme_pubclean(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=13),
        axis.text.x  = element_text(size=13),
        axis.title.y = element_text(size=13),
        axis.text.y  = element_text(vjust=0.5, size=13) ) +
  theme(legend.text = element_text(size = 10), legend.title = element_text(size=10), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=15.5)) +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  theme(legend.position = "bottom")

xaxislevelsF1 <- c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog", "genic","antisense","fusion","intergenic","genic_intron");
xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")
subc.levels=c("alternative_3end",'alternative_3end5end', "alternative_5end","reference_match", "3prime_fragment","internal_fragment", "5prime_fragment","combination_of_known_junctions", "combination_of_known_splicesites", "intron_retention","no_combination_of_known_junctions", "mono-exon_by_intron_retention", "at_least_one_novel_splicesite", "mono-exon", "multi-exon")
subc.labels=c("Alternative 3'end", "Alternative 3'5'end", "Alterantive 5'end", "Reference match", "3' fragment", "Internal fragment", "5' fragment", "Comb. of annot. junctions", "Comb. of annot. splice sites", "Intron retention", "Not comb. of annot. junctions", "Mono-exon by intron ret.", "At least 1 annot. don./accept.", "Mono-exon", "Multi-exon")
coding.levels=c("coding", "non_coding")
coding.labels=c("Coding", "Non coding")

rename_SC_labels <- function(x){
  x$structural_category = factor(x$structural_category,
                                 labels = xaxislabelsF1,
                                 levels = xaxislevelsF1,
                                 ordered=TRUE)
  x$subcategory = factor(x$subcategory,
                         labels = subc.labels,
                         levels = subc.levels,
                         ordered=TRUE)
  x$coding = factor(x$coding,
                    labels = coding.labels,
                    levels = coding.levels,
                    ordered=TRUE)
  x
}

rename_SC_labels_sirvs <- function(x){
  x$structural_category = factor(x$structural_category,
                                 labels = xaxislabelsF1,
                                 levels = xaxislevelsF1,
                                 ordered=TRUE)
   x
}

raw_classif <- read.csv("SourceData_figures/HIS/WTC11_cDNA.HIS_classification.txt", header=T, sep="\t")

raw_classif <- rename_SC_labels(raw_classif)

legendLabelF1 <- levels(as.factor(raw_classif$coding))

fig2A <- ggplot(data=raw_classif, aes(x=structural_category)) +
  geom_bar(aes(y = (..count..)/sum(..count..)*100, alpha=coding, fill=structural_category), color="black", width=0.7, size=0.3) +
  geom_text(aes(y = (..count..)/sum(..count..)*100, label = (..count..)), stat = "count", vjust = -0.25)  +
  scale_x_discrete(drop=FALSE) +
  scale_alpha_manual(values=c(1,0.3),
                     name = "Coding prediction",
                     labels = legendLabelF1)+
  xlab("") +
  ylab("Transcripts, %") +
  pub_theme +
  theme(axis.text.x = element_text(angle = 45, hjust=0.8, vjust = 1.2)) +
  scale_fill_manual(values = cat.palette, guide='none') +
  ggtitle("Isoform Distribution Across Structural Categories\n" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12)) +
  scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))

ggsave(filename = "main_figures/figure2A.pdf", plot = fig2A, width=6, height=6)


fig2B <- ggplot(data=raw_classif, aes(x=within_CAGE_peak, y=ratio_TSS, fill=within_CAGE_peak)) +
  geom_half_violin(alpha=0.8, adjust=1,
                   width=1, side="r" )+
  geom_half_boxplot(width = 1,
                    outlier.shape = NA, side="l" ) +
  stat_compare_means(method = "wilcox.test",
                     label = "p.signif",
                     label.y = 2.5, label.x=1.5) +
  scale_x_discrete(drop=TRUE) +
  pub_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5)) +
  scale_fill_conesa(palette="warm", name="", reverse = T) +
  labs(x="\nDetected TSS \n\nwithin a CAGE peak",  y="TSS ratio") +
  scale_y_continuous(expand = expansion(mult=c(0,0.1)),
                     trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)),
                     limits = c(0.1,500)) +
  theme(legend.position="none") +
  theme(axis.title.x = element_text(lineheight=.4, size=15))

ggsave(filename = "main_figures/figure2B.pdf", plot = fig2B, width=3, height=6)

## Main Figure 2C


FSM_ISM_classif <- raw_classif  %>% filter(structural_category %in% c("FSM", "ISM")) %>% 
  mutate(less200bp_TSS = ifelse(abs(diff_to_gene_TSS)<=50,
                                "Known",
                                "Novel"),
         low_ratio_TSS = ifelse(ratio_TSS<1.5,
                                "TSS ratio < 1.5",
                                "TSS ratio > 1.5"))

fig2c.1 <- ggplot(data=FSM_ISM_classif, aes(x=ratio_TSS, fill=structural_category)) +
  geom_density(alpha=0.7) +
  geom_vline(xintercept=1.5, linetype="dashed", color = "red")+
  ggh4x::facet_grid2( .~ within_CAGE_peak , scales="free",  independent = "y",
                      labeller = as_labeller(c(`TRUE`="Supported by CAGE-seq", 
                                               `FALSE`="Not supported by CAGE-seq")) ) +
  scale_x_continuous(expand = expansion(mult=c(0,0.1)),
                     trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  pub_theme +
  theme(axis.text.x = element_text(angle=0)) +
  scale_fill_manual(values = cat.palette) +
  #scale_fill_conesa(palette="main", name="TSS annotated") +
  labs(y="Density", x="TSS ratio") +
  theme(legend.position = "none") +
  guides(fill=guide_legend(title="")) 

fig2c.2 <- ggplot(FSM_ISM_classif, aes(y=structural_category, x=..count..,
                                        alpha=less200bp_TSS, fill=structural_category)) +
  geom_bar(color="black", position="stack") +
  facet_wrap(.~ within_CAGE_peak, labeller = as_labeller(c(`TRUE`="Supported by CAGE-seq", 
                                                           `FALSE`="Not supported by CAGE-seq")))+
  pub_theme +
  theme_pubclean(flip=T)+
  labs(x="Total number of isoforms",
       y="") +
  scale_fill_manual(values = cat.palette) +
  scale_alpha_manual(values = c(1,0.4)) +
  #scale_fill_conesa("warm", reverse=T, name="") +
  scale_x_continuous(label = unit_format(unit = "K", scale = 0.001, accuracy = 1), expand = expansion(mult=c(0,0.2))) + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size=12), legend.title = element_text(size=12)) +
  guides(fill=guide_legend(title=""),
         alpha=guide_legend(title = "TSS status"))


fig2c <- fig2c.1 / fig2c.2 + 
  plot_layout(heights = c(4, 1), ncol = 1)

ggsave(filename = "main_figures/figure2C.pdf", plot = fig2c, width = 6.5, height = 5.5)

### Figure 2D: ML importance

ML_importance_df <- read.csv("SourceData_figures/HIS/ML_filter/classifier_variable-importance_table.txt", sep="\t", header = F)

fig42d <- ggplot(ML_importance_df, aes(y=V2, x=reorder(V1, -V2), fill=""))+
  geom_bar(stat = "identity", position = "dodge", color="black")+
  pub_theme +
  scale_fill_conesa(palette="main", name="", reverse = T ) +
  scale_y_continuous(expand = expansion(mult=c(0,0.01)))+
  theme(axis.text.x = element_text(size=13, angle=30, hjust=1)) +
  xlab("") + ylab("") +
  labs(title = "ML filter",
       subtitle = "Variable importance") +
  theme(legend.position="none")


ggsave(filename = "main_figures/figure2D.pdf", plot = fig42d, width=7, height=5)

## Figure 2e: Value distribution for ML relevant variants across structural categories

class_ML_HIS <- read.csv("SourceData_figures/HIS/ML_filter/WTC11_cDNA.HIS_MLresult_classification.txt", sep="\t",
                         header = T, as.is=T)
class_ML_HIS <- rename_SC_labels(class_ML_HIS)
reduced_classML <- class_ML_HIS %>%
  select(isoform, chrom, structural_category, subcategory,
         ratio_TSS, diff_to_gene_TSS, FL, min_cov, iso_exp,
         diff_to_gene_TTS, filter_result) %>% filter(!str_detect(chrom,"ERCC|SIRV"))

pivoted_reduced_class <- reduced_classML %>% pivot_longer(cols = c("ratio_TSS", "diff_to_gene_TSS","FL",
                                                                   "min_cov","iso_exp","diff_to_gene_TTS"))

pivoted_reduced_class <- reduced_classML %>% pivot_longer(cols = c("ratio_TSS", "diff_to_gene_TSS","FL",
                                                                   "min_cov"))


pivoted_reduced_class$name <- pivoted_reduced_class$name %>% 
  factor(levels=c("ratio_TSS", "diff_to_gene_TSS","FL",
                  "min_cov","iso_exp","diff_to_gene_TTS"),
         labels = c("TSS ratio", "Distance to\nknown TSS (bp)",
                    "Num. FL reads", "Minimum Short-Read\ncoverage",
                    "Isoform expression (TPM)", "Distance to\nknown TTS (bp)"))

ylims_adjusted <- pivoted_reduced_class %>% na.omit() %>% 
  group_by(name) %>%
  summarise(Q1 = quantile(value, 1/4), Q3 = quantile(value, 3/4)) %>%
  ungroup() %>%
  #get lowest Q1 and highest Q3
  summarise(lowQ1 = 0, highQ3 = log10(max(Q3)+1))

fig2e <- ggplot(pivoted_reduced_class %>% filter(structural_category %in% c("FSM", "ISM","NIC","NNC")) ,
            aes(y=log10(value+1), alpha=filter_result, fill=structural_category,
                x=structural_category))+
  geom_boxplot(width=1,outlier.shape = NA) +
  scale_fill_manual(values = cat.palette) +
  scale_alpha_manual(values = c(0.1, 0.8))+
  facet_wrap(name~., scales = "free", nrow = 2, drop = T ) +
  pub_theme +
  theme(strip.text = element_text(size=10),
        axis.title.x = element_blank())+
  guides(alpha = guide_legend(override.aes = list(fill = "#5a5b5c"), 
                              title=""),
         fill = guide_legend(title=""))

fig2e <- fig2e + coord_cartesian(ylim = as.numeric(ylims_adjusted)*2 ) 

ggsave(filename = "main_figures/figure2E.pdf", plot = fig2e, width=7, height=4)

## Figure 2F

rescue_names_HI1 <- c("SourceData_figures/HIS/rules_filter/rescue/WTC11_cDNA.HIS_rules_rescue_classification.txt",
                      "SourceData_figures/HIS/ML_filter/rescue/WTC11_cDNA.HIS_rescue_classification.txt")
names <- c("rules_rescue", "ML_rescue")

HI_1 <- list()
HI_1[["raw"]]<-raw_classif
HI_1[["rules_filter"]] <-read.csv("SourceData_figures/HIS/rules_filter/WTC11_cDNA.HIS_RulesFilter_result_classification.txt", sep="\t", header = T) %>%
               rename_SC_labels() %>% 
               filter(filter_result=="Isoform")
HI_1[["ML_filter"]]<-read.csv("SourceData_figures/HIS/ML_filter/WTC11_cDNA.HIS_MLresult_classification.txt", sep="\t", header = T) %>%
               rename_SC_labels() %>% 
               filter(filter_result=="Isoform")

for (j in 1:2){
  HI_1[[names[j]]] <- read.csv(rescue_names_HI1[j], sep="\t", header = T) %>% rename_SC_labels()
}

get_num_isof_knGene <- function(x){
  y <- x %>% filter(structural_category %in% c("FSM", "ISM", "NIC", "NNC")) %>% select(isoform, associated_gene)
  r <- data.frame("Type info"=c("Num Isoforms", "Num Known Genes"),
                  "Value"=c(y$isoform %>% length() , 
                            y$associated_gene %>% unique() %>%  length())
  )
  r
}

absolut_numbers_HI_1 <- map(HI_1, get_num_isof_knGene) %>% bind_rows(.id = "Stage")

all_num <- bind_rows(list("High Input Sample"=absolut_numbers_HI_1),
                     .id="Input")

all_num$Stage_pipeline <- all_num$Stage %>% factor( levels=c("raw", "ML_filter", "ML_rescue", "rules_filter", "rules_rescue"),
                                                    labels=c("Pre-filter", "Filter", "Rescue", "Filter", "Rescue"))

all_num$Type_pipeline <- all_num$Stage %>% factor( levels=c("raw", "ML_filter", "ML_rescue", "rules_filter", "rules_rescue"),
                                                   labels=c("ML", "ML", "ML", "Rules", "Rules"))

add_initial_rules <- all_num %>% filter(Stage=="raw")
add_initial_rules$Type_pipeline = "Rules"

all_num <- rbind(all_num, add_initial_rules)

all_num2 <- all_num %>% mutate(newInput=ifelse(Stage=="raw",
                                               "Starting transcriptome",
                                               Input)) %>% select(-Input) %>% unique()

all_num2$Type_pipeline <- all_num2$Type_pipeline %>% as.character()
all_num2[which(all_num2$Stage=="raw"), "Type_pipeline"]="IsoSeq3"
all_num2 <- all_num2 %>% unique()

all_num2$Type_pipeline <- all_num2$Type_pipeline %>% factor(levels=c("IsoSeq3", "ML", "Rules"),
                                                            labels=c("IsoSeq3", "IsoSeq + SQ3-ML",
                                                                     "IsoSeq + SQ3-Rules"))

fig2f <- ggplot(all_num2,
                    aes(x=Stage_pipeline, y=Value, fill=Type_pipeline)) +
  geom_bar(stat="identity", position = position_dodge2(preserve = "single"), color="black")+
  ggh4x::facet_grid2( .~Type.info , scales="free"  , independent = "y") +
  pub_theme +
  scale_fill_met_d("Egypt") +
  theme(legend.position="bottom") +
  theme(axis.text.x=element_text(angle=0, vjust=1, hjust=0.5),
        strip.text = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_y_continuous(labels = unit_format(unit = "K", scale = 0.001, accuracy = 1),  
                     expand = expansion(mult=c(0,0.1)), 
                     position = "left") +
  guides(fill=guide_legend(title=""))


ggsave(filename = "main_figures/figure2F.pdf", plot = fig2f, width=8, height=4)

### SIRVs plot


HI_1_SIRVs <- list()
HI_1_SIRVs[["raw"]]<- read.csv("SourceData_figures/HIS/SIRVs/SIRV_mod_class.raw.extended.txt", sep="\t", header = F) 
colnames(HI_1_SIRVs[["raw"]]) <- c("isoform","chrom", "structural_category", "associated_transcript",	"diff_to_TSS",	"diff_to_TTS",
                             "all_canonical",	"min_sample_cov",	"min_cov",	"FL",	"iso_exp",	"gene_exp",	"ratio_exp",	"ratio_TSS")

HI_1_SIRVs[["raw"]] <- HI_1_SIRVs[["raw"]] %>% rename_SC_labels_sirvs()
colnames(HI_1_SIRVs[["raw"]]) <- paste0(colnames(HI_1_SIRVs[["raw"]]), ".Mod")

HI_1_SIRVs[["rules_filter"]]<- read.csv("SourceData_figures/HIS/SIRVs/SIRV_mod_class.rules.extended.txt", sep="\t", header = F)
colnames(HI_1_SIRVs[["rules_filter"]]) <- c("isoform","chrom", "structural_category", "associated_transcript",	"diff_to_TSS",	"diff_to_TTS",
                             "all_canonical",	"min_sample_cov",	"min_cov",	"FL",	"iso_exp",	"gene_exp",	"ratio_exp",	"ratio_TSS")
HI_1_SIRVs[["rules_filter"]] <- HI_1_SIRVs[["rules_filter"]] %>% rename_SC_labels_sirvs()
colnames(HI_1_SIRVs[["rules_filter"]]) <- paste0(colnames(HI_1_SIRVs[["rules_filter"]]), ".Mod")

HI_1_SIRVs[["ML_filter"]]<- read.csv("SourceData_figures/HIS/SIRVs/SIRV_mod_class.ML.extended.txt", sep="\t", header = F) 
colnames(HI_1_SIRVs[["ML_filter"]]) <- c("isoform","chrom", "structural_category", "associated_transcript",	"diff_to_TSS",	"diff_to_TTS",
                                            "all_canonical",	"min_sample_cov",	"min_cov",	"FL",	"iso_exp",	"gene_exp",	"ratio_exp",	"ratio_TSS")
HI_1_SIRVs[["ML_filter"]] <- HI_1_SIRVs[["ML_filter"]] %>% rename_SC_labels_sirvs()
colnames(HI_1_SIRVs[["ML_filter"]]) <- paste0(colnames(HI_1_SIRVs[["ML_filter"]]), ".Mod")

HI_1_SIRVs[["rules_rescue"]]<- read.csv("SourceData_figures/HIS/SIRVs/SIRV_mod_class.rules_rescue.extended.txt", sep="\t", header = F)
colnames(HI_1_SIRVs[["rules_rescue"]]) <- c("isoform","chrom", "structural_category", "associated_transcript",	"diff_to_TSS",	"diff_to_TTS",
                                            "all_canonical",	"min_sample_cov",	"min_cov",	"FL",	"iso_exp",	"gene_exp",	"ratio_exp",	"ratio_TSS")
HI_1_SIRVs[["rules_rescue"]] <- HI_1_SIRVs[["rules_rescue"]] %>% rename_SC_labels_sirvs()
colnames(HI_1_SIRVs[["rules_rescue"]]) <- paste0(colnames(HI_1_SIRVs[["rules_rescue"]]), ".Mod")

HI_1_SIRVs[["ML_rescue"]]<- read.csv("SourceData_figures/HIS/SIRVs/SIRV_mod_class.ML_rescue.extended.txt", sep="\t", header = F) 
colnames(HI_1_SIRVs[["ML_rescue"]]) <- c("isoform","chrom", "structural_category", "associated_transcript",	"diff_to_TSS",	"diff_to_TTS",
                                            "all_canonical",	"min_sample_cov",	"min_cov",	"FL",	"iso_exp",	"gene_exp",	"ratio_exp",	"ratio_TSS")
HI_1_SIRVs[["ML_rescue"]] <- HI_1_SIRVs[["ML_rescue"]] %>% rename_SC_labels_sirvs()
colnames(HI_1_SIRVs[["ML_rescue"]]) <- paste0(colnames(HI_1_SIRVs[["ML_rescue"]]), ".Mod")


sirv_list=read.table("SourceData_figures/HIS/SIRVs/SIRVs_ids.txt", header = F)$V1
sirv_list <- sirv_list[grep(pattern = "00",x = sirv_list , invert = T )]

HI_1_SIRVs_mod <- map2(HI_1, HI_1_SIRVs, merge, by.x="isoform", by.y="isoform.Mod")

get_hits <- function(x){
  x$hit <- apply(x,1,function(y){
    ifelse(y["structural_category"]=="FSM" & 
             abs(as.numeric(y["diff_to_TSS"]))<=50 & abs(as.numeric(y["diff_to_TTS"]))<=50 & 
             y["associated_transcript.Mod"] %in% sirv_list,
           "known_TP",
           ifelse(y["structural_category"]=="FSM" & 
                    abs(as.numeric(y["diff_to_TSS"]))<=50 & abs(as.numeric(y["diff_to_TTS"]))<=50 &
                    ! y["associated_transcript.Mod"] %in% sirv_list,
                  "novel_TP",
                  ifelse(! y["associated_transcript.Mod"] %in% c(sirv_list, "novel"),
                         "overAnnot_FP",
                         ifelse(y["structural_category"] %in% c("FSM", "ISM") &
                                  (abs(as.numeric(y["diff_to_TSS"]))>50 | abs(as.numeric(y["diff_to_TTS"]))>50),
                                "shortened_TP",
                                "FP"))))
  })
  x
}

HI_1_SIRVs_hit <- map(HI_1_SIRVs_mod, get_hits)


get_performance <- function(x){
  # get rough numbers
  num_sirvs <- length(sirv_list)
  sirvs_called <- length(x$isoform)
  matched_sirvs_called <- x %>% filter(structural_category %in% c("FSM", "ISM")) %>% select(isoform)
  TP_called <- x %>% filter(hit %in% c("known_TP","novel_TP")) %>% 
    select(associated_transcript) %>% unique() 
  novel_TP_called <-  x %>% filter(hit=="novel_TP") %>% 
    select(associated_transcript) %>% 
    unique() 
  known_TP_called <-  x %>% filter(hit=="known_TP") %>% 
    select(associated_transcript) %>% 
    unique()
  PTP_called <- x %>% filter(hit=="shortened_TP") %>% 
    select(associated_transcript)
  FP_called <- x %>% filter(hit=="FP") %>% select(isoform)
  FN_called <- setdiff(sirv_list, TP_called$associated_transcript)
  
  #get metrics
  Sensitivity=length(TP_called$associated_transcript)*100/num_sirvs %>% round(digits = 2)
  PDR=length(unique(c(TP_called$associated_transcript, PTP_called$associated_transcript)))*100/num_sirvs %>% round(digits = 2)
  Precision=length(TP_called$associated_transcript)*100/sirvs_called %>% round(digits = 2)
  redPrecision=(length(TP_called$associated_transcript)+length(PTP_called$associated_transcript))*100/sirvs_called %>% round(digits = 2)
  FDiscR=(length(FP_called$isoform) + length(PTP_called$associated_transcript))*100/sirvs_called %>% round(digits = 2)
  FDetR=(length(FP_called$isoform))*100/sirvs_called %>% round(digits = 2)
  redPrecision=(100-FDetR) %>% round(digits = 2)
  Redundancy = length(matched_sirvs_called$isoform) / length(unique(c(TP_called$associated_transcript, PTP_called$associated_transcript)))
  
  Fscore <- (2*Sensitivity*Precision)/(Sensitivity+Precision)
  
  # using mod columns, we can evaluate the novel discovery rate and missleading effect of overannotation
  novel_SIRVs <- x %>% filter(! structural_category %in% c("FSM", "ISM"))
  NDR <- length(novel_SIRVs$isoform)*100/sirvs_called %>% round(digits = 2)
  
  overannot_called <- x %>% filter(hit=="overAnnot_FP") 
  unique_OvFP <- overannot_called %>% select(associated_transcript) %>% unique()
  
  ODR <- length(overannot_called$isoform)*100/sirvs_called %>% round(digits = 2)
  
  result <- data.frame("Total detected"=sirvs_called,
                       "TP"=length(TP_called$associated_transcript), "PTP"=length(PTP_called$associated_transcript),
                       "novel_TP"=length(novel_TP_called$associated_transcript),
                       "known_TP"=length(known_TP_called$associated_transcript),
                       "FP"=length(FP_called$isoform), "FN"=length(FN_called),
                       "ovFP"=length(overannot_called$isoform),
                       "unique_ov_FP"=length(unique_OvFP$associated_transcript),
                       "Sensitivity"=Sensitivity, "Precision"=Precision, "RedPrecision"=redPrecision, "Fscore"=Fscore,
                       "PDR"=PDR, "FDiscR"=FDiscR, "FDetR"=FDetR, "NDR"=NDR, "ODR"=ODR, "Redundancy"=Redundancy)
  result <- apply(result,2, round, 2) %>% as.data.frame() %>% t() %>% as.data.frame()
  return(result)
}

performance_HI_1 <- map(HI_1_SIRVs_hit, get_performance)

performance_HI_1_df <- bind_rows(performance_HI_1, .id = "step")

rownames(performance_HI_1_df)<- 1:length(performance_HI_1_df$step)

performance_all_df <- bind_rows(list("High Input Sample"=performance_HI_1_df),
                                .id = "Input")

get_stage<-function(x){
  if (x["step"] %in% c("rules_filter","ML_filter")){
    return("Filter")
  }
  if (x["step"] %in% c("rules_rescue","ML_rescue")){
    return("Rescue")
  }
  if (x["step"]=="raw"){
    return("Pre-filter")
  }
}

get_type_1<-function(x){
  if (x["step"] %in% c("rules_filter","rules_rescue")){
    return("Rules")
  }
  if (x["step"] %in% c("ML_filter","ML_rescue")){
    return("ML")
  }
  if (x["step"]=="raw"){
    return("ML")
  }
}

get_type <- function(x){
  x$Type <- apply(x,1,get_type_1)
  to_add <- x %>% filter(step=="raw")
  to_add$Type <- "Rules"
  y <- rbind(x,to_add)
  return(y)
}
performance_all_df$Stage <- apply(performance_all_df, 1, get_stage) %>% factor(levels=c("Pre-filter","Filter","Rescue"))
performance_all_df <- get_type(performance_all_df)


## plot 

padf <- performance_all_df %>% filter(Input=="High Input Sample") %>% 
pivot_longer(cols = c(Sensitivity,Precision,Fscore,FDiscR, NDR, ODR)) %>%
  mutate(metric=factor(.$name , levels = c("Sensitivity","Precision","Fscore","FDiscR","NDR","ODR"),
                       labels = c("Sensitivity","Precision","F-score","FDR","NDR","ODR")))

p_sirvs <- ggplot(padf, aes(x=Stage, y=value)) +
  geom_segment( aes(x=Stage, xend=Stage, y=0, yend=value, color=metric), size=1.3) +
  geom_point( size=2.5, aes(color=metric ))  +
  geom_line(aes(group=metric, color=metric), linetype=2, size=0.7) +
  facet_grid( Type ~ metric, scales = "free", space = "free", switch = "y"  ) +
  pub_theme +
  scale_color_conesa(palette = "main") +
  labs(x="", y="") +
  theme(legend.position="none",
        axis.text.x = element_text(angle=90, vjust=0.5, hjust = 1)) +
  scale_y_continuous(expand=expansion(mult=c(0,0.1)),limits = c(0, 100),
                     labels = scales::percent_format(scale = 1), position="right")

ggsave(filename = "main_figures/figure2G.pdf", plot = p_sirvs, width=8, height=4.5)






