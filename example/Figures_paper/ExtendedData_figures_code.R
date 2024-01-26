library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(RColorConesa)
library(scales)
library(patchwork)
library(ggpubr)

if (!file.exists("ExtendedData_figures")) {
  dir.create("ExtendedData_figures")
}

cat.palette = c( "FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                 "NNC"="#EE6A50", "GenicGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                 "Intergenic" = "darksalmon", "GenicIntron"="#41B6C4")

cat.palette2 = c( "FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                 "NNC"="#EE6A50", "Genic\nGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                 "Intergenic" = "darksalmon", "Genic\nIntron"="#41B6C4")

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

raw_classif <- read.csv("SourceData_figures/HIS/WTC11_cDNA.HIS_classification.txt", header=T, sep="\t")

raw_classif <- rename_SC_labels(raw_classif)

legendLabelF1 <- levels(as.factor(raw_classif$coding))

## read HIR classification to know CAGE bs refTSS overlap

## Extended Data Figure 2
# VennDiagram Overlap TSS 
hir_info <- read.csv("SourceData_figures/HIR/WTC11_cDNA.HIR_classification.txt", header=T, sep="\t") %>% 
  select(isoform, within_CAGE_peak, within_polyA_site)

colnames(hir_info) <- c("isoform", "within_refTSS_peak","within_atlas_polyA_peak")

df_class2 <- merge(raw_classif, hir_info, by="isoform")

a <- table(df_class2$within_CAGE_peak, df_class2$within_refTSS_peak, as.numeric(df_class2$ratio_TSS)>=1.5) %>% as.data.frame()
colnames(a) <- c("CAGE", "refTSS", "ratio_TSS", "Freq")

library(venneuler)
MyVenn <- venneuler(c("A"=a %>% filter(CAGE==T, refTSS==F, ratio_TSS==F) %>% select(Freq) %>% as.numeric(),
                      "B"=a %>% filter(CAGE==F, refTSS==T, ratio_TSS==F) %>% select(Freq) %>% as.numeric(),
                      "C"=a %>% filter(CAGE==F, refTSS==F, ratio_TSS==T) %>% select(Freq) %>% as.numeric(),
                      "A&B"=a %>% filter(CAGE==T, refTSS==T, ratio_TSS==F) %>% select(Freq) %>% as.numeric(), 
                      "A&C"=a %>% filter(CAGE==T, refTSS==F, ratio_TSS==T) %>% select(Freq) %>% as.numeric(),
                      "B&C"=a %>% filter(CAGE==F, refTSS==T, ratio_TSS==T) %>% select(Freq) %>% as.numeric(),
                      "A&B&C"= a %>% filter(CAGE==T, refTSS==T, ratio_TSS==T) %>% select(Freq) %>% as.numeric()
                      ))

MyVenn$labels <- ""

cage_num <- a %>% filter(CAGE==T) %>% select(Freq) %>% sum()
cage_num <- as.numeric(cage_num)*0.001 
cage_num <- cage_num %>% round(digits=0)

refTSS_num <- a %>% filter(refTSS==T) %>% select(Freq) %>% sum()
refTSS_num <- as.numeric(refTSS_num)*0.001 
refTSS_num <- refTSS_num %>% round(digits=0)

ratioTSS_num <- a %>% filter(ratio_TSS==T) %>% select(Freq) %>% sum()
ratioTSS_num <- as.numeric(ratioTSS_num)*0.001 
ratioTSS_num <- ratioTSS_num %>% round(digits=0)



all_num=a %>% filter(CAGE==T, refTSS==T, ratio_TSS==T) %>% select(Freq) %>% as.numeric()
all_num <- as.numeric(all_num)*0.001 
all_num <- all_num %>% round(digits=1)

pdf("ExtendedData_figures/ExtendedData_Figure2.pdf")
plot(MyVenn)

text(MyVenn$centers[1, 1] + 0.21, MyVenn$centers[1, 2] + 0.24,
     paste0("CAGE-seq\nn=",cage_num,"K"), cex = 2)
text(MyVenn$centers[2, 1] + 0.21, MyVenn$centers[2, 2] - 0.24,
     paste0("refTSS\nn=", refTSS_num, "K") , cex = 2)
text(MyVenn$centers[3, 1] - 0.24, MyVenn$centers[3, 2] + 0.24, 
     paste0("TSS ratio\nn=",ratioTSS_num,"K"), cex = 2)
text(0.5, 0.5, paste0(all_num,"K\ncalled by all"), cex = 2)

dev.off()

## Extended Data Figure 3
# Upset Plot Overlap TTS 

validated_TTS_list <- list()

validated_TTS_list[["Quant-seq"]] <- df_class2 %>% filter(within_polyA_site==TRUE)
validated_TTS_list[["Quant-seq"]] <- validated_TTS_list[["Quant-seq"]]$isoform

validated_TTS_list[["polyASite"]] <- df_class2 %>% filter(within_atlas_polyA_peak==TRUE)
validated_TTS_list[["polyASite"]] <- validated_TTS_list[["polyASite"]]$isoform

validated_TTS_list[["polyA motif"]] <- df_class2 %>% filter(polyA_motif_found==TRUE)
validated_TTS_list[["polyA motif"]] <- validated_TTS_list[["polyA motif"]]$isoform

library(UpSetR)
library(ComplexHeatmap)

m_TTS = list_to_matrix(validated_TTS_list) %>% as.data.frame()
m_TTS$isoform <- rownames(m_TTS)
m_TTS<- merge(m_TTS, df_class2, by="isoform")


upsetPlot_TTS <- upset(m_TTS,
                       query.legend = "bottom",
                       order.by = "freq",
                       text.scale = c(c(2, 1.5, 1.5, 1, 1.5, 1.5)),
                       #c(1.3, 1.3, 1, 1, 2, 0.75),
                       queries = list(
                         list(
                           query = elements,
                           params=list("structural_category", c("FSM", "ISM", "NIC","NNC")),
                           color = "#6BAED6",
                           active = T,
                           query.name = "FSM",
                           query.size=12
                         ),
                         list(
                           query = elements,
                           params=list("structural_category", c("ISM", "NIC","NNC")),
                           color = "#FC8D59",
                           active = T,
                           query.name = "ISM",
                           query.size=12
                         ),
                         list(
                           query = elements,
                           params=list("structural_category", c("NIC","NNC")),
                           color = "#78C679",
                           active = T,
                           query.name = "NIC",
                           query.size=12
                         ),
                         list(
                           query = elements,
                           params=list("structural_category", "NNC"),
                           color = "#EE6A50",
                           active = T,
                           query.name = "NNC",
                           query.size=12
                         )
                       )
)

pdf("ExtendedData_figures/ExtendedData_Figure3.pdf")
upsetPlot_TTS
dev.off()


## Extended Data Figure 4
# QuantSeq vs polyA motif distance

quantseq_status <- c("TRUE"="Supported by Quant-seq",
                     "FALSE"="Not supported by Quant-seq")

supp_fig4 <- ggplot(data=raw_classif %>% filter(polyA_motif_found==T), aes( x=polyA_dist, fill=structural_category)) +
  geom_histogram(binwidth = 2) +
  geom_vline(xintercept = -17, linetype="dashed") +
  scale_fill_manual(values = cat.palette, limits=force, name="") +
  facet_wrap( ~ within_polyA_site, scales = "free", 
             labeller = labeller(structural_category=label_value,
                                 within_polyA_site=quantseq_status)) +
  pub_theme +
  theme(legend.position = "bottom") 

ggsave(filename = "ExtendedData_figures/ExtendedData_Figure4.pdf", plot = supp_fig4, width=8, height=5)

## intraapriming vs ployA motif distance
raw_classif_intr <- raw_classif %>%  mutate(intrapriming=ifelse(perc_A_downstream_TTS>=60,
                                                                TRUE,
                                                                FALSE))

freq <- table(raw_classif_intr$structural_category, raw_classif_intr$intrapriming)

raw_classif_intr <- raw_classif_intr %>%
  group_by(structural_category, intrapriming) %>%
  mutate(count_per_group = n())

raw_classif_intr$count_per_group <- as.factor(paste0("n=", raw_classif_intr$count_per_group))

ord <- raw_classif_intr %>% select(structural_category, count_per_group, intrapriming) %>%
  unique() %>% as.data.frame() %>% arrange(structural_category)

raw_classif_intr$count_per_group <- raw_classif_intr$count_per_group %>% factor(
  levels=ord$count_per_group, labels=ord$count_per_group
)
  
supp_fig5 <-ggplot(raw_classif_intr,
       aes(y = polyA_dist, x = (count_per_group), fill = structural_category)) +
  geom_violin() + 
  geom_boxplot(width = 0.3, color = "lightgrey", outlier.shape = NA, show.legend = F) +
  facet_wrap(.~ intrapriming, scales = "free", labeller = label_both) +
  pub_theme +
  scale_fill_manual(values = cat.palette2, limits = force, name = "") +
  geom_hline(yintercept = -17, color = "black", linetype = "dashed") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1) )+
  labs(x = "Count per Group",
       y = "Distance between motif\nand polyA site")


ggsave(filename = "ExtendedData_figures/ExtendedData_Figure5.pdf", plot = supp_fig5, width=8, height=5)

# plot tiles

classif <- read.csv("SourceData_figures/HIS/ML_filter/rescue/WTC11_cDNA.HIS_rescue_classification.txt", sep="\t", header = T) %>% rename_SC_labels()
classif$structural_category <- classif$structural_category %>%  as.character()
classif[grep(classif$isoform,pattern = "ENST"), "structural_category"] <- "Reference"
classif[grep(classif$isoform,pattern = "DQ"), "structural_category"] <- "Reference"
classif[grep(classif$isoform,pattern = "SIRV"), "structural_category"] <- "Reference"
classif[grep(classif$isoform,pattern = "EF"), "structural_category"] <- "Reference"

classif$structural_category <- classif$structural_category %>% as.factor()

rescue_table <- read.csv("SourceData_figures/HIS/ML_filter/rescue/WTC11_cDNA.HIS_rescue_table.tsv", sep = "\t", header=T)

rescue_table[which(rescue_table$sam_flag==4),"exclusion_reason"]<- "Unmapped"
rescue_table[which(rescue_table$rescue_result=="rescued_automatic"),"exclusion_reason"]<- "Automatic rescue"
rescue_table[which(rescue_table$rescue_result=="rescued_mapping"),"exclusion_reason"]<- "Mapping rescue"


rescue_table_SC <- merge(rescue_table, classif[, c("isoform", "structural_category")],
                         by.x="mapping_hit", by.y="isoform", all.x = T)
rescue_table_SC[which(rescue_table_SC$sam_flag==4), "structural_category"]<- "Unmapped"
rescue_table_SC[is.na(rescue_table_SC$structural_category), "structural_category"]<- "Reference"

colnames(rescue_table_SC) <- c("mapping_hit", "rescue_candidate", "sam_flag", "POS_MLprob", "structural_category_candidate",
                               "Rescue_result","exclusion_reason","best_match_for_candidate", "best_match_id","structural_category_hit")

### make frequency tables

first_table <- table(rescue_table_SC$structural_category_candidate, rescue_table_SC$structural_category_hit) %>% as.data.frame()
first_table <- first_table[first_table$Freq!=0, ]

colnames(first_table) <- c("Source", "Target", "Freq")

first_table$Source <- first_table$Source %>%  factor(levels = c("full-splice_match",
                                                                "incomplete-splice_match",
                                                                "novel_in_catalog",
                                                                "novel_not_in_catalog"),
                                                     labels=c(paste0("FSM discarded\nn=",
                                                                     rescue_table_SC %>% filter(structural_category_candidate =="full-splice_match") %>%
                                                                       select(rescue_candidate) %>% 
                                                                       unlist() %>% unique() %>% length()),
                                                              paste0("ISM discarded\nn=",
                                                                     rescue_table_SC %>% filter(structural_category_candidate =="incomplete-splice_match") %>%
                                                                       select(rescue_candidate) %>% 
                                                                       unlist() %>% unique() %>% length()),
                                                              paste0("NIC discarded\nn=",
                                                                     rescue_table_SC %>% filter(structural_category_candidate =="novel_in_catalog") %>%
                                                                       select(rescue_candidate) %>% 
                                                                       unlist() %>% unique() %>% length()),
                                                              paste0("NNC discarded\nn=",
                                                                     rescue_table_SC %>% filter(structural_category_candidate =="novel_not_in_catalog") %>%
                                                                       select(rescue_candidate) %>% 
                                                                       unlist() %>% unique() %>% length()))
)

first_table$Target <- first_table$Target %>%  factor(levels = c("Reference",
                                                                "FSM",
                                                                "ISM",
                                                                "NIC",
                                                                "NNC",
                                                                "Fusion",
                                                                "Unmapped"),
                                                     labels=c(paste0("Reference hit\nn=",
                                                                     rescue_table_SC %>% filter(structural_category_hit =="Reference") %>%
                                                                       select(mapping_hit) %>% 
                                                                       unlist() %>% unique() %>% length()),
                                                              paste0("FSM hit\nn=",
                                                                     rescue_table_SC %>% filter(structural_category_hit =="FSM") %>%
                                                                       select(mapping_hit) %>% 
                                                                       unlist() %>% unique() %>% length()),
                                                              paste0("ISM hit\nn=",
                                                                     rescue_table_SC %>% filter(structural_category_hit =="ISM") %>%
                                                                       select(mapping_hit) %>% 
                                                                       unlist() %>% unique() %>% length()),
                                                              paste0("NIC hit\nn=",
                                                                     rescue_table_SC %>% filter(structural_category_hit =="NIC") %>%
                                                                       select(mapping_hit) %>% 
                                                                       unlist() %>% unique() %>% length()),
                                                              paste0("NNC hit\nn=",
                                                                     rescue_table_SC %>% filter(structural_category_hit =="NNC") %>%
                                                                       select(mapping_hit) %>% 
                                                                       unlist() %>% unique() %>% length()),
                                                              paste0("Fusion hit\nn=",
                                                                     rescue_table_SC %>% filter(structural_category_hit =="Fusion") %>%
                                                                       select(mapping_hit) %>% 
                                                                       unlist() %>% unique() %>% length()),
                                                              paste0("Unmapped\nn=",
                                                                     rescue_table_SC %>% filter(structural_category_hit =="Unmapped") %>%
                                                                       select(rescue_candidate) %>% 
                                                                       unlist() %>% unique() %>% length()))
                                                     
)

first_table$Freq <- first_table$Freq %>% as.numeric()
SC_Source <- str_split(first_table$Source, pattern=" ") %>% unlist()
row_odd <- seq_len(length(SC_Source)) %% 2
SC_Source <- SC_Source[row_odd==1]
first_table$SC_Source <- SC_Source

Nsource <- data.frame("SC_Source"=c("FSM","ISM","NIC","NNC"),
                      "Num_source"=c(rescue_table_SC %>% filter(structural_category_candidate =="full-splice_match") %>%
                                       select(rescue_candidate) %>% 
                                       unlist() %>% unique() %>% length(),
                                     rescue_table_SC %>% filter(structural_category_candidate =="incomplete-splice_match") %>%
                                       select(rescue_candidate) %>% 
                                       unlist() %>% unique() %>% length(),
                                     rescue_table_SC %>% filter(structural_category_candidate =="incomplete-splice_match") %>%
                                       select(rescue_candidate) %>% 
                                       unlist() %>% unique() %>% length(),
                                     rescue_table_SC %>% filter(structural_category_candidate =="novel_not_in_catalog") %>%
                                       select(rescue_candidate) %>% 
                                       unlist() %>% unique() %>% length()
                      ))

first_table <- merge(first_table, Nsource, by="SC_Source") 
first_table$hits_per_candidate <- first_table$Freq/first_table$Num_source 
first_table$hits_per_candidate  <- first_table$hits_per_candidate %>% round(digits = 2)

first_table$label <- paste0(first_table$Freq, "\n(",first_table$hits_per_candidate, ")")
plot.tiles <- ggplot(first_table, aes(y=Source, x=Target, fill=log10(Freq) )) +
  geom_tile()+
  geom_text(aes(label=label)) +
  scale_fill_viridis(discrete = F, name="log10(Num. hits)") +
  pub_theme +
  theme(legend.justification=c(1,1), legend.direction = "horizontal",
        legend.position = c(1,0.2)) +
  scale_x_discrete(position = "top") +
  labs(x="Unique Targets",
       y="Rescue Candidates")

ggsave(filename = "ExtendedData_figures/ExtendedData_Figure6.pdf", plot = plot.tiles, width=9, height=7)





