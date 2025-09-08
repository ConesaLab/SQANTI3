#####################################
#####       TUSCO Report       ######
#####################################
### Author: Tianyuan Liu
### Last Modified: 11/07/2024 by tianyuan.liu@csic.es

#********************** Taking arguments from python script
args <- commandArgs(trailingOnly = TRUE)
class.file <- args[1]
tusco.file <- args[2]
transcript_gtf_file <- args[3]
utilities.path <- args[4]
# optional genome assembly (e.g., hg38, mm10)
genome.assembly <- ifelse(length(args) >= 5, args[5], "hg38")

# Define file paths

# class.file <- '/media/tian/ubuntu/SQANTI_TUSCO/WTC11_cdna_ont/WTC11_cdna_ont_classification.txt'
# tusco.file <- '/media/tian/ubuntu/GitHub/SQANTI3/utilities/report_qc/tusco_human.gtf'
# transcript_gtf_file <- '/media/tian/ubuntu/SQANTI_TUSCO/HUMAN/WTC11_cdna_ont.gtf'
# utilities.path <- '/media/tian/ubuntu/GitHub/SQANTI3/utilities'

report.prefix <- strsplit(class.file, "_classification.txt")[[1]][1]
output_directory <- dirname(class.file)
output_name <- basename(report.prefix)
html.report.file <- paste0(output_name, "_TUSCO_report.html")

# Load necessary libraries
message("Loading necessary libraries...")
suppressPackageStartupMessages({
  library(ggplot2)
  library(rmarkdown)
  library(tidyr)
  library(dplyr)
  library(readr)
  library(stringr)
  library(rtracklayer)
  library(GenomicRanges)
  library(Gviz)
})

# Helper function to read TSV files with error handling
read_tsv_safe <- function(file_path, col_names = TRUE, ...) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  tryCatch(
    {
      data <- read_tsv(file_path, col_names = col_names, ...)
      message("Successfully read file: ", file_path)
      return(data)
    },
    error = function(e) {
      stop("Error reading file ", file_path, ": ", e$message)
    }
  )
}

# Read input files
classification_data <- read_tsv_safe(class.file)

## Read TUSCO reference (supports TSV or GTF)
if (grepl("\\.tsv$", tusco.file, ignore.case = TRUE) || grepl("\\.txt$", tusco.file, ignore.case = TRUE)) {
  message("Reading TUSCO TSV file...")
  tusco_df <- tryCatch({
    readr::read_delim(
      tusco.file,
      delim = "\t",
      col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
      col_types = readr::cols(.default = "c"),
      trim_ws = TRUE
    )
  }, error = function(e) {
    # fallback parser if needed
    readr::read_table2(
      tusco.file,
      col_names = c("ensembl", "transcript", "gene_name", "gene_id_num", "refseq", "prot_refseq"),
      col_types = readr::cols(.default = "c")
    )
  })
  if (!all(c("ensembl","refseq","gene_name") %in% colnames(tusco_df))) {
    stop("Tusco TSV must contain columns: ensembl, refseq, gene_name")
  }
  annotation_data <- tusco_df %>%
    dplyr::select(ensembl, refseq, gene_name) %>%
    dplyr::distinct()
  # For downstream plotting code, keep a unified object name
  tusco_gtf_df <- tusco_df
} else {
  message("Reading TUSCO GTF file...")
  tusco_gtf <- rtracklayer::import(tusco.file)
  tusco_gtf_df <- as.data.frame(tusco_gtf)

  # Derive the expected columns if missing
  if (!"ensembl" %in% names(tusco_gtf_df)) {
    if ("gene_id" %in% names(tusco_gtf_df)) {
      tusco_gtf_df$ensembl <- tusco_gtf_df$gene_id
    } else {
      tusco_gtf_df$ensembl <- NA_character_
    }
  }
  if (!"gene_name" %in% names(tusco_gtf_df)) {
    tusco_gtf_df$gene_name <- NA_character_
  }
  if (!"refseq" %in% names(tusco_gtf_df)) {
    tusco_gtf_df$refseq <- NA_character_
  }

  # Extract gene-level information. Some GTFs might lack explicit 'gene' rows,
  # so derive unique identifiers from any feature rows available.
  annotation_data <- tusco_gtf_df %>%
    dplyr::select(ensembl = ensembl, refseq = refseq, gene_name = gene_name) %>%
    dplyr::filter(!is.na(ensembl) | !is.na(refseq) | !is.na(gene_name)) %>%
    dplyr::distinct()
}

## initial rTUSCO (may be recomputed after selecting id type)
rTUSCO <- nrow(annotation_data)

# Read transcript GTF file
message("Reading transcript GTF file...")
transcript_gtf <- rtracklayer::import(transcript_gtf_file)
transcript_gtf_df <- as.data.frame(transcript_gtf)
message("transcript_gtf_df columns: ", paste(colnames(transcript_gtf_df), collapse = ", "))
message("tusco_gtf_df columns: ", paste(colnames(tusco_gtf_df), collapse = ", "))

message("Defining regex patterns for ID classification...")
patterns <- list(
  ensembl   = "^(ENSG|ENSMUSG)\\d{11}(\\.\\d+)?$",
  refseq    = "^(NM_|NR_|NP_)\\d{6,}$",
  gene_name = "^[A-Z0-9]+$"
)

## Defer summary printing until id_summary is computed

message("Cleaning the classification data...")
# Clean the classification data (split fusions, strip versions), then classify id_type
message("Cleaning the classification data...")
classification_data_cleaned <- classification_data %>%
  filter(structural_category != "fusion") %>%
  bind_rows(
    classification_data %>%
      filter(structural_category == "fusion") %>%
      separate_rows(associated_gene, sep = "_")
  ) %>%
  mutate(
    associated_gene = str_remove(associated_gene, "\\.\\d+$")
  ) %>%
  distinct(isoform, associated_gene, .keep_all = TRUE) %>%
  arrange(isoform) %>%
  mutate(
    id_type = case_when(
      str_detect(associated_gene, patterns$ensembl)   ~ "ensembl",
      str_detect(associated_gene, patterns$refseq)    ~ "refseq",
      str_detect(associated_gene, patterns$gene_name) ~ "gene_name",
      TRUE ~ "unknown"
    )
  )

## Generate and display summary of id_type counts from cleaned data
id_summary <- dplyr::count(classification_data_cleaned, id_type, sort = TRUE)
message("Summary of id_type counts:")
print(id_summary)

# Identify and display the id_type with the highest count
if (nrow(id_summary) > 0) {
  top_id_row <- id_summary %>% dplyr::slice_max(n, n = 1)
  top_id_type <- top_id_row$id_type
  top_id_n <- top_id_row$n
  message("The id_type with the highest count is: ", top_id_type, " with ", top_id_n, " entries.")
} else {
  message("No id_type classifications were made.")
}

# Normalize IDs: strip version suffixes from relevant fields
if (top_id_type == "ensembl" && "gene_id" %in% names(transcript_gtf_df)) {
  # Directly modify the 'gene_id' column
  transcript_gtf_df$gene_id <- sub("\\..*", "", transcript_gtf_df$gene_id)
}

# Also normalize annotation IDs to remove version suffixes
if ("ensembl" %in% names(annotation_data)) {
  annotation_data$ensembl <- sub("\\..*", "", as.character(annotation_data$ensembl))
}
if ("refseq" %in% names(annotation_data)) {
  annotation_data$refseq <- sub("\\..*", "", as.character(annotation_data$refseq))
}

# Normalize TUSCO GTF identifiers used for matching/plotting
if ("ensembl" %in% names(tusco_gtf_df)) {
  tusco_gtf_df$ensembl <- sub("\\..*", "", as.character(tusco_gtf_df$ensembl))
}
if ("refseq" %in% names(tusco_gtf_df)) {
  tusco_gtf_df$refseq <- sub("\\..*", "", as.character(tusco_gtf_df$refseq))
}

## associated_transcript version strip (if present) after cleaning
if ("associated_transcript" %in% names(classification_data_cleaned)) {
  classification_data_cleaned$associated_transcript <- stringr::str_remove(classification_data_cleaned$associated_transcript, "\\\\.\\\\d+$")
}

# Metrics and definitions for evaluation against TUSCO
message("Defining TUSCO-related transcript mappings...")
# TUSCO_transcripts: Transcripts mapping to TUSCO genes (by selected ID type)
annotation_ids <- annotation_data[[top_id_type]]
annotation_ids <- unique(annotation_ids[!is.na(annotation_ids)])
TUSCO_transcripts <- classification_data_cleaned %>%
  filter(associated_gene %in% annotation_ids)

# TUSCO_RM: True Positives
# Rule: (1) reference_match OR (2) FSM mono-exon with both ends within 50bp
TUSCO_RM <- TUSCO_transcripts %>%
  filter(
    subcategory == "reference_match" |
    (structural_category == "full-splice_match" & !is.na(ref_exons) & ref_exons == 1 &
       !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
       abs(diff_to_TSS) <= 50 & abs(diff_to_TTS) <= 50)
  )

# Define True Positives (TP): TUSCO transcripts identified as Reference Match (RM)
TP <- TUSCO_RM
TP_TUSCO <- unique(TUSCO_RM$associated_gene)

# Define Partial True Positives (PTP): TUSCO transcripts identified as FSM or ISM but not RM
PTP <- TUSCO_transcripts %>%
  filter(structural_category %in% c("full-splice_match", "incomplete-splice_match") & !associated_gene %in% TP$associated_gene)

# Define False Negatives (FN): TUSCO genes with no SQANTI3 transcript in any category
FN <- annotation_data %>%
  filter(!(!!sym(top_id_type) %in% TUSCO_transcripts$associated_gene))

# Define False Positives (FP): Transcripts in NIC, NNC, genic, or fusion categories within TUSCO_transcripts
FP <- TUSCO_transcripts %>%
  filter(structural_category %in% c("novel_in_catalog", "novel_not_in_catalog", "genic", "fusion"))

fsm_ism_count <- TUSCO_transcripts %>%
  filter(structural_category %in% c("full-splice_match", "incomplete-splice_match")) %>%
  nrow()

# Remove only unneeded columns; keep schema even if empty
TP <- TP %>% select(-associated_transcript, -id_type)
PTP <- PTP %>% select(-associated_transcript, -id_type)
FP <- FP %>% select(-associated_transcript, -id_type)

# Calculate metrics
message("Calculating evaluation metrics...")

# Initialize variables to store metric values
non_redundant_sensitivity <- NA
non_redundant_precision <- NA
redundant_precision <- NA
positive_detection_rate <- NA
false_discovery_rate <- NA
false_detection_rate <- NA
redundancy <- NA

# Calculate Non-redundant Precision: TP / TUSCO_transcripts
if (nrow(TUSCO_transcripts) > 0) {
  non_redundant_precision <- nrow(TP) / nrow(TUSCO_transcripts)
} else {
  warning("TUSCO_transcripts has zero rows. Non-redundant Precision set to NA.")
}

# Calculate Precision: (TP + PTP) / TUSCO_transcripts
message("Calculating redundant Precision...")
if (fsm_ism_count > 0) {
  redundant_precision <- (nrow(TP) + nrow(PTP)) /  nrow(TUSCO_transcripts)
  message("Precision calculated successfully.")
} else {
  redundant_precision <- NA
  warning("FSM + ISM count is zero. Redundant Precision set to NA.")
}

# Calculate Sensitivity: TP_TUSCO/rTUSCO
message("Calculating Sensitivity...")
non_redundant_sensitivity <- length(TP_TUSCO)/rTUSCO

# Calculate Positive Detection Rate: unique(TP + PTP) / rTUSCO
unique_detected_genes <- length(unique(c(TP$associated_gene, PTP$associated_gene)))
if (rTUSCO > 0) {
  positive_detection_rate <- unique_detected_genes / rTUSCO
} else {
  warning("annotation_data has zero rows. Positive Detection Rate set to NA.")
}

# Calculate False Discovery Rate: (TUSCO_transcripts - TUSCO_RM) / TUSCO_transcripts
if (nrow(TUSCO_transcripts) > 0) {
  false_discovery_rate <- (nrow(TUSCO_transcripts) - nrow(TUSCO_RM)) / nrow(TUSCO_transcripts)
} else {
  warning("TUSCO_transcripts has zero rows. False Discovery Rate set to NA.")
}

# Calculate False Detection Rate: FP / TUSCO_transcripts
if (nrow(TUSCO_transcripts) > 0) {
  false_detection_rate <- nrow(FP) / nrow(TUSCO_transcripts)
} else {
  warning("TUSCO_transcripts has zero rows. False Detection Rate set to NA.")
}
# Calculate Redundancy: (FSM + ISM) / unique(TP + PTP)
unique_tp_ptp_genes <- length(unique(c(TP$associated_gene, PTP$associated_gene)))
if (unique_tp_ptp_genes > 0) {
  redundancy <- fsm_ism_count / unique_tp_ptp_genes
} else {
  warning("No unique TP + PTP genes found. Redundancy set to NA.")
}

# Compile metrics into a data frame for better readability
metrics <- tibble(
  Metric = c(
    "Sensitivity",
    "Non-redundant Precision",
    "Redundant Precision",
    "Positive Detection Rate",
    "False Discovery Rate",
    "False Detection Rate",
    "Redundancy"
  ),
  Value = c(
    non_redundant_sensitivity,
    non_redundant_precision,
    redundant_precision,
    positive_detection_rate,
    false_discovery_rate,
    false_detection_rate,
    redundancy
  )
) %>%
  mutate(
    Value = round(Value * 100, 2),
    Value = paste0(Value, "%")
  )

# Display the metrics
message("Evaluation Metrics:")
print(metrics)

# Prepare data for plotting
message("Preparing data for plotting...")
# Combine TP, PTP, and FP datasets
combined_data <- bind_rows(
  TP %>% mutate(category = "TP"),
  PTP %>% mutate(category = "PTP"),
  FP %>% mutate(category = "FP")
)

# Extract unique associated_genes for each category
associated_genes_list <- combined_data %>%
  dplyr::select(associated_gene, category) %>%
  dplyr::distinct()

# Include FN genes for plotting (reference-only tracks if no sample transcripts)
fn_ids <- tryCatch({
  if (nrow(FN) > 0) as.character(FN[[top_id_type]]) else character(0)
}, error = function(e) character(0))
genes_to_plot <- unique(c(associated_genes_list$associated_gene, fn_ids))
genes_to_plot <- genes_to_plot[!is.na(genes_to_plot) & genes_to_plot != ""]
message("IGV genes to plot: ", length(genes_to_plot))

## NOTE: params will be created after all derived objects are computed
####################################################################
# TUSCO: SQANTI structural category
####################################################################
# Filter TUSCO_only to keep transcripts that match TUSCO genes
# TUSCO_only should already be defined as TUSCO_transcripts.
TUSCO_only <- TUSCO_transcripts

# Rename subcategories to match your palettes
TUSCO_only <- TUSCO_only %>%
  mutate(subcategory = case_when(
    subcategory == "reference_match" ~ "Reference match",
    subcategory == "alternative_3end" ~ "Alternative 3'end",
    subcategory == "alternative_5end" ~ "Alternative 5'end",
    subcategory == "alternative_3end5end" ~ "Alternative 3'5'end",
    TRUE ~ subcategory
  ))

# Define big categories:
# a (TP): FSM + Reference match = "RM"
# b (PTP): FSM not RM (Alternative ends) + ISM
# c (FP): NIC, NNC, Genic Intron, Genic Genomic, Antisense, Fusion, Intergenic
# d (FN): Missing genes

TUSCO_only <- TUSCO_only %>%
  mutate(
    mono_exon_close50 = structural_category == "full-splice_match" & !is.na(ref_exons) & ref_exons == 1 &
                        !is.na(diff_to_TSS) & !is.na(diff_to_TTS) &
                        abs(diff_to_TSS) <= 50 & abs(diff_to_TTS) <= 50,
    big_category = case_when(
      structural_category == "full-splice_match" & (subcategory == "Reference match" | mono_exon_close50) ~ "TP",
      (structural_category == "full-splice_match" & subcategory %in% c("Alternative 3'end", "Alternative 5'end", "Alternative 3'5'end")) ~ "PTP",
      structural_category == "incomplete-splice_match" ~ "PTP",
      structural_category %in% c("novel_in_catalog","novel_not_in_catalog","genic_intron","genic","antisense","fusion","intergenic") ~ "FP",
      TRUE ~ NA_character_
    )
  )

# final_label for plotting
TUSCO_only <- TUSCO_only %>%
  mutate(final_label = case_when(
    big_category == "TP" ~ "RM",
    big_category == "PTP" & subcategory == "Alternative 3'end" ~ "Alternative 3'end",
    big_category == "PTP" & subcategory == "Alternative 5'end" ~ "Alternative 5'end",
    big_category == "PTP" & subcategory == "Alternative 3'5'end" ~ "Alternative 3'5'end",
    big_category == "PTP" & structural_category == "incomplete-splice_match" ~ "ISM",
    big_category == "FP" & structural_category == "novel_in_catalog" ~ "NIC",
    big_category == "FP" & structural_category == "novel_not_in_catalog" ~ "NNC",
    big_category == "FP" & structural_category == "genic_intron" ~ "Genic Intron",
    big_category == "FP" & structural_category == "genic" ~ "Genic Genomic",
    big_category == "FP" & structural_category == "antisense" ~ "Antisense",
    big_category == "FP" & structural_category == "fusion" ~ "Fusion",
    big_category == "FP" & structural_category == "intergenic" ~ "Intergenic",
    TRUE ~ NA_character_
  ))

# Create data frame for FN
missing_df <- data.frame(
  final_label = "Missing",
  big_category = "FN",
  count = nrow(FN)
)

# Summarize TUSCO_only counts
plot_data <- TUSCO_only %>%
  filter(!is.na(final_label)) %>%
  group_by(big_category, final_label) %>%
  summarise(count = n(), .groups="drop")

# Add missing data
plot_data <- bind_rows(plot_data, missing_df)

# Calculate percentages
plot_data <- plot_data %>%
  mutate(percentage = count/sum(count)*100)

# Palettes
cat.palette <- c("FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                 "NNC"="#EE6A50", "Genic Genomic"="#969696", "Antisense"="#66C2A4", 
                 "Fusion"="goldenrod1","Intergenic"="darksalmon", "Genic Intron"="#41B6C4")

subcat.palette <- c("Alternative 3'end"='#02314d',
                    "Alternative 3'5'end"='#0e5a87',
                    "Alternative 5'end"='#7ccdfc',
                    "Reference match"='#c4e1f2',
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

# Map final_label to colors:
# a: RM (use "Reference match" color)
# b: Alternative ends (from subcat.palette), ISM (cat.palette["ISM"])
# c: NIC, NNC, Genic Intron, Genic Genomic, Antisense, Fusion, Intergenic (from cat.palette)
# d: Missing (grey60)

final_colors <- c(
  "RM" = subcat.palette[["Reference match"]],
  "Alternative 3'end" = subcat.palette[["Alternative 3'end"]],
  "Alternative 5'end" = subcat.palette[["Alternative 5'end"]],
  "Alternative 3'5'end" = subcat.palette[["Alternative 3'5'end"]],
  "ISM" = cat.palette[["ISM"]],
  "NIC" = cat.palette[["NIC"]],
  "NNC" = cat.palette[["NNC"]],
  "Genic Intron" = cat.palette[["Genic Intron"]],
  "Genic Genomic" = cat.palette[["Genic Genomic"]],
  "Antisense" = cat.palette[["Antisense"]],
  "Fusion" = cat.palette[["Fusion"]],
  "Intergenic" = cat.palette[["Intergenic"]],
  "Missing" = "grey60"
)


# Set factor levels
plot_data$big_category <- factor(plot_data$big_category, levels=c("TP","PTP","FP","FN"))
plot_data$final_label <- factor(plot_data$final_label,
                                levels = c("RM",
                                           "Alternative 3'end","Alternative 5'end","Alternative 3'5'end","ISM",
                                           "NIC","NNC","Genic Intron","Genic Genomic","Antisense","Fusion","Intergenic",
                                           "Missing"))

# Theme
mytheme <- theme_classic(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4),
        axis.title.x = element_text(size=13),
        axis.text.x  = element_text(size=12, angle=45, hjust=1),
        axis.title.y = element_text(size=13),
        axis.text.y  = element_text(vjust=0.5, size=12),
        legend.position="none",
        plot.title = element_text(lineheight=.4, size=15, hjust = 0.5))

# Plot
p_tusco_complex <- ggplot(plot_data, aes(x = final_label, y = percentage, fill = final_label)) +
  geom_bar(stat="identity", color="black", size=0.3, width=0.7) +
  xlab(NULL) +
  facet_grid(~big_category, scales="free_x", space="free") +
  scale_fill_manual(values = final_colors) +
  scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
  mytheme

# Pass data to RMarkdown (after all objects are defined)
params <- list(
  metrics = metrics,
  tusco_gtf_df = tusco_gtf_df,
  transcript_gtf_df = transcript_gtf_df,
  associated_genes_list = associated_genes_list,
  p_tusco_complex = p_tusco_complex,
  # pass objects to avoid global leakage
  TP = TP,
  PTP = PTP,
  FN = FN,
  FP = FP,
  non_redundant_sensitivity = non_redundant_sensitivity,
  non_redundant_precision = non_redundant_precision,
  redundant_precision = redundant_precision,
  positive_detection_rate = positive_detection_rate,
  false_discovery_rate = false_discovery_rate,
  redundancy = redundancy,
  output_directory = output_directory,
  genome_assembly = genome.assembly
)
########################################################
# Generate IGV plots for each gene
########################################################

message("Generating IGV plots for each gene...")

# Create a directory to save the plots
plots_dir <- file.path(output_directory, "igv_plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}
# Improve headless PNG rendering (e.g., on servers/CI)
options(ucscChromosomeNames = FALSE)
if (tolower(getOption("bitmapType", default = "")) != "cairo") {
  try({ options(bitmapType = "cairo") }, silent = TRUE)
}

# Map transcript_id to associated_gene for coloring/labels (keep original gene_id too)
if ("transcript_id" %in% names(transcript_gtf_df)) {
  transcript_gtf_df$gene_id_original <- transcript_gtf_df$gene_id
  tidx <- match(transcript_gtf_df$transcript_id, classification_data_cleaned$isoform)
  map_ok <- !is.na(tidx)
  message("Mapped transcript_id to classification isoform: ", sum(map_ok), "/", length(map_ok), " matched")
  # Only remap gene_id if we actually matched sample isoforms; skip for reference GTF
  if (any(map_ok)) {
    transcript_gtf_df$gene_id[map_ok] <- classification_data_cleaned$associated_gene[tidx[map_ok]]
  }
  # strip version from gene_id for matching
  if ("gene_id" %in% names(transcript_gtf_df)) {
    transcript_gtf_df$gene_id <- sub("\\.\\d+$", "", as.character(transcript_gtf_df$gene_id))
  }
}
plots_created <- 0L
if (length(genes_to_plot) == 0) {
  message("No genes selected for IGV plotting (genes_to_plot is empty).")
}
msg_preview <- paste(utils::head(genes_to_plot, 15), collapse = ", ")
message("IGV plotting preview (up to 15 genes): ", msg_preview)
message("plots_dir: ", plots_dir)
message("capabilities(cairo)=", paste0(isTRUE(capabilities("cairo"))), ", getOption(bitmapType)=", as.character(getOption("bitmapType")))
for (gene in genes_to_plot) {
  # Fetch TUSCO gene annotations (exons) for the current gene, if available
  tusco_has_coords <- all(c("seqnames","start","end","strand") %in% names(tusco_gtf_df))
  tusco_gene <- tryCatch({
    if (tusco_has_coords) {
      .df <- tusco_gtf_df
      if ("type" %in% names(.df)) {
        .df <- dplyr::filter(.df, type == "exon")
      }
      # Build a safe match vector that ignores NAs in matching columns
      match_vec <- rep(FALSE, nrow(.df))
      if ("gene_name" %in% names(.df)) match_vec <- match_vec | (!is.na(.df$gene_name) & as.character(.df$gene_name) == gene)
      if ("ensembl" %in% names(.df))   match_vec <- match_vec | (!is.na(.df$ensembl)   & as.character(.df$ensembl)   == gene)
      if ("refseq" %in% names(.df))    match_vec <- match_vec | (!is.na(.df$refseq)    & as.character(.df$refseq)    == gene)
      .df <- .df[match_vec %in% TRUE, , drop = FALSE]
      # Drop rows with missing coordinates proactively
      dplyr::filter(.df, !is.na(seqnames) & !is.na(start) & !is.na(end) & !is.na(strand))
    } else {
      # Derive TUSCO gene coordinates from the user-provided reference GTF
      .ref <- transcript_gtf_df
      if ("type" %in% names(.ref)) {
        .ref <- dplyr::filter(.ref, type == "exon")
      }
      # Build match sets from TUSCO TSV annotations (prefer gene_name, then ensembl)
      tusco_gene_names <- annotation_data$gene_name[!is.na(annotation_data$gene_name)]
      tusco_ensembl_ids <- annotation_data$ensembl[!is.na(annotation_data$ensembl)]
      tusco_ensembl_ids <- sub("\\.\\d+$", "", as.character(tusco_ensembl_ids))
      match_vec <- rep(FALSE, nrow(.ref))
      if ("gene_name" %in% names(.ref) && length(tusco_gene_names) > 0) {
        match_vec <- match_vec | (!is.na(.ref$gene_name) & .ref$gene_name %in% tusco_gene_names)
      }
      if ("gene_id" %in% names(.ref) && length(tusco_ensembl_ids) > 0) {
        gid <- sub("\\.\\d+$", "", as.character(.ref$gene_id))
        match_vec <- match_vec | (!is.na(gid) & gid %in% tusco_ensembl_ids)
      }
      .ref <- .ref[match_vec %in% TRUE, , drop = FALSE]
      dplyr::filter(.ref, !is.na(seqnames) & !is.na(start) & !is.na(end) & !is.na(strand))
    }
  }, error = function(e) data.frame())
  
  # Fetch transcript exons associated with this gene
  gene_transcripts <- tryCatch({
    .tdf <- transcript_gtf_df
    if ("type" %in% names(.tdf)) {
      .tdf <- dplyr::filter(.tdf, type == "exon")
    }
    match_vec <- rep(FALSE, nrow(.tdf))
    if ("gene_id" %in% names(.tdf)) match_vec <- match_vec | (!is.na(.tdf$gene_id)   & sub("\\\\.\\\\d+$", "", as.character(.tdf$gene_id)) == gene)
    if ("gene_name" %in% names(.tdf)) match_vec <- match_vec | (!is.na(.tdf$gene_name) & as.character(.tdf$gene_name) == gene)
    .tdf <- .tdf[match_vec %in% TRUE, , drop = FALSE]
    dplyr::filter(.tdf, !is.na(seqnames) & !is.na(start) & !is.na(end) & !is.na(strand))
  }, error = function(e) transcript_gtf_df[0, , drop = FALSE])
  
  message("Gene ", gene, ": tusco_exons=", nrow(tusco_gene), ", sample_exons=", nrow(gene_transcripts))
  if (nrow(gene_transcripts) > 0) {
    message("Sample seqnames: ", paste(unique(as.character(gene_transcripts$seqnames)), collapse=","))
  }
  if (nrow(tusco_gene) > 0) {
    message("TUSCO seqnames: ", paste(unique(as.character(tusco_gene$seqnames)), collapse=","))
  }

  message(paste("Creating plot for gene:", gene))
  
  tryCatch({
    message("Starting plot generation for gene: ", gene)
    if (nrow(gene_transcripts) == 0 && nrow(tusco_gene) == 0) {
      message("No transcript or reference features found for ", gene, "; skipping plot.")
      next
    }
    
    # Merge gene_transcripts with classification_data_cleaned using transcript_id and isoform
    # Merge classification info; continue gracefully if join fails
    gene_transcripts <- tryCatch({
      labels_df <- classification_data_cleaned %>%
        dplyr::select(isoform, structural_category) %>%
        dplyr::mutate(read_name = isoform)
      suppressWarnings(
        gene_transcripts %>%
          left_join(labels_df, by = c("transcript_id" = "isoform"))
      )
    }, error = function(e) {
      message("left_join failed: ", e$message, "; proceeding without join.")
      gene_transcripts
    })
    
    # Map structural_category to abbreviations
    category_abbr <- c(
      "full-splice_match" = "FSM",
      "incomplete-splice_match" = "ISM",
      "novel_in_catalog" = "NIC",
      "novel_not_in_catalog" = "NNC",
      "genic" = "Genic\nGenomic",
      "genic_intron" = "Genic\nIntron",
      "intergenic" = "Intergenic",
      "antisense" = "Antisense",
      "fusion" = "Fusion"
    )
    gene_transcripts$structural_category_abbr <- category_abbr[gene_transcripts$structural_category]
    
    # Map abbreviations to colors
    cat.palette <- c(
      "FSM" = "#6BAED6",
      "ISM" = "#FC8D59",
      "NIC" = "#78C679",
      "NNC" = "#EE6A50",
      "Genic\nGenomic" = "#969696",
      "Antisense" = "#66C2A4",
      "Fusion" = "goldenrod1",
      "Intergenic" = "darksalmon",
      "Genic\nIntron" = "#41B6C4"
    )
    gene_transcripts$color <- cat.palette[gene_transcripts$structural_category_abbr]

    # Build a mapping from transcript_id -> desired track label (prefer SQANTI isoform/read name)
    # Fallback to transcript_id if isoform is missing
    name_map <- tryCatch({
      track_names_df <- gene_transcripts %>% dplyr::select(transcript_id, read_name) %>% dplyr::distinct()
      # Coalesce missing isoform to transcript_id
      track_names_df$track_label <- ifelse(is.na(track_names_df$read_name) | track_names_df$read_name == "",
                                           as.character(track_names_df$transcript_id),
                                           as.character(track_names_df$read_name))
      setNames(track_names_df$track_label, track_names_df$transcript_id)
    }, error = function(e) {
      message("Failed to build name_map: ", e$message)
      NULL
    })
    
    sample_transcripts_gr_unlisted <- NULL
    if (nrow(gene_transcripts) > 0) {
      message("Creating GRangesList for sample transcripts...")
      sample_transcripts_gr <- makeGRangesListFromDataFrame(
        gene_transcripts,
        keep.extra.columns = TRUE,
        seqnames.field = "seqnames",
        start.field = "start",
        end.field = "end",
        strand.field = "strand",
        split.field = "transcript_id"
      )
      sample_transcripts_gr_unlisted <- unlist(sample_transcripts_gr, use.names = TRUE)
      mcols(sample_transcripts_gr_unlisted)$transcript <- rep(
        names(sample_transcripts_gr),
        elementNROWS(sample_transcripts_gr)
      )
      # Assign colors to mcols
      transcript_colors <- gene_transcripts %>%
        dplyr::select(transcript_id, color) %>%
        dplyr::distinct()
      color_map <- setNames(transcript_colors$color, transcript_colors$transcript_id)
      mcols(sample_transcripts_gr_unlisted)$fill <- color_map[mcols(sample_transcripts_gr_unlisted)$transcript]
      mcols(sample_transcripts_gr_unlisted)$col <- color_map[mcols(sample_transcripts_gr_unlisted)$transcript]
    }
    
    tusco_transcripts_gr_unlisted <- NULL
    if (tusco_has_coords && nrow(tusco_gene) > 0) {
      # Create GRangesList for TUSCO reference transcripts
      message("Creating GRangesList for TUSCO reference transcripts...")
      split_field <- if ("transcript_id" %in% names(tusco_gene)) "transcript_id" else "ensembl"
      tusco_transcripts_gr <- makeGRangesListFromDataFrame(
        tusco_gene,
        keep.extra.columns = TRUE,
        seqnames.field = "seqnames",
        start.field = "start",
        end.field = "end",
        strand.field = "strand",
        split.field = split_field
      )
      tusco_transcripts_gr_unlisted <- unlist(tusco_transcripts_gr, use.names = TRUE)
      mcols(tusco_transcripts_gr_unlisted)$transcript <- rep(
        names(tusco_transcripts_gr),
        elementNROWS(tusco_transcripts_gr)
      )
      # Set fill and col to 'darkgrey' for TUSCO
      mcols(tusco_transcripts_gr_unlisted)$fill <- 'darkgrey'
      mcols(tusco_transcripts_gr_unlisted)$col <- 'darkgrey'
    }
    
    # Define the plotting range
    message("Defining the plotting range...")
    if (!is.null(tusco_transcripts_gr_unlisted) && nrow(tusco_gene) > 0) {
      gene_range_start <- min(c(gene_transcripts$start, tusco_gene$start))
      gene_range_end <- max(c(gene_transcripts$end, tusco_gene$end))
      gene_seqnames <- unique(c(gene_transcripts$seqnames, tusco_gene$seqnames))
      gene_strand <- unique(c(gene_transcripts$strand, tusco_gene$strand))
    } else {
      gene_range_start <- min(gene_transcripts$start, na.rm = TRUE)
      gene_range_end <- max(gene_transcripts$end, na.rm = TRUE)
      gene_seqnames <- unique(gene_transcripts$seqnames)
      gene_strand <- unique(gene_transcripts$strand)
    }
    gene_range <- GRanges(
      seqnames = gene_seqnames[1],
      ranges = IRanges(start = gene_range_start, end = gene_range_end),
      strand = gene_strand[1]
    )
    message("gene_range: chr=", as.character(seqnames(gene_range)[1]), ", start=", gene_range_start, ", end=", gene_range_end, ", strand=", as.character(strand(gene_range)[1]))
    
    genome_axis <- GenomeAxisTrack()
    message("Genome axis created.")
    
    # Generate GeneRegionTrack for sample transcripts (if any)
    sample_tracks <- list()
    if (!is.null(sample_transcripts_gr_unlisted) && length(sample_transcripts_gr_unlisted) > 0) {
      message("Creating GeneRegionTrack for sample transcripts...")
      message("Sample transcripts to track: ", length(unique(mcols(sample_transcripts_gr_unlisted)$transcript)))
      sample_tracks <- lapply(unique(mcols(sample_transcripts_gr_unlisted)$transcript), function(transcript) {
        # Subset GRanges for the current transcript
        current_transcript_gr <- sample_transcripts_gr_unlisted[mcols(sample_transcripts_gr_unlisted)$transcript == transcript]
        
        # Fetch the color for the current transcript
        current_color <- unique(mcols(current_transcript_gr)$fill)
        # Determine the display label for the track: prefer SQANTI isoform (read name)
        label <- as.character(transcript)
        if (!is.null(name_map) && as.character(transcript) %in% names(name_map)) {
          nm <- name_map[[as.character(transcript)]]
          if (!is.na(nm) && nzchar(nm)) label <- nm
        }
        
        # Create a GeneRegionTrack for the current transcript
        GeneRegionTrack(
          current_transcript_gr,
          genome = genome.assembly,
          chromosome = as.character(seqnames(gene_range)[1]),
          name = label,
          background.title = current_color,
          fill = current_color,
          col = current_color
        )
      })
    } else {
      message("No sample transcripts for ", gene, "; plotting TUSCO reference only if available.")
    }
    
    # Create GeneRegionTrack for TUSCO reference
    tusco_track <- NULL
    if (!is.null(tusco_transcripts_gr_unlisted)) {
      message("Creating GeneRegionTrack for TUSCO reference...")
      message("TUSCO reference transcripts to track: ", length(unique(mcols(tusco_transcripts_gr_unlisted)$transcript)))
      tusco_track <- GeneRegionTrack(
        tusco_transcripts_gr_unlisted,
        genome = genome.assembly,
        chromosome = as.character(seqnames(gene_range)[1]),
        name = paste0(gene, " TUSCO Reference"),
        background.title = "darkgrey",
        fill = "darkgrey",
        col = "darkgrey"
      )
    }
    # Combine all transcript tracks into a single list of tracks
    sample_tracks <- if (!is.null(sample_transcripts_gr_unlisted)) unlist(sample_tracks, recursive = FALSE) else list()
    
    # Combine all tracks for plotting
    if (is.null(tusco_track)) {
      all_tracks <- c(genome_axis, sample_tracks)
    } else {
      all_tracks <- c(genome_axis, tusco_track, sample_tracks)
    }
    if (length(all_tracks) <= 1) {
      message("No tracks to plot for gene: ", gene)
      next
    }
    message("All tracks combined for plotting.")
    
    # Plot and save to file
    plot_file <- file.path(plots_dir, paste0(gene, ".png"))
    message("Saving plot to file: ", plot_file)
    # Prefer Cairo device when available to avoid X11 issues on headless systems
    open_png <- function(file) {
      ok <- FALSE
      if (isTRUE(capabilities("cairo"))) {
        try({ png(file, width = 1000, height = 600, type = "cairo"); ok <- TRUE }, silent = TRUE)
      }
      if (!ok) {
        png(file, width = 1000, height = 600)
      }
    }
    message("Opening PNG device (cairo=", isTRUE(capabilities("cairo")), ", bitmapType=", as.character(getOption("bitmapType")), ")")
    open_png(plot_file)
    plotTracks(
      all_tracks,
      from = gene_range_start,
      to = gene_range_end,
      chromosome = as.character(seqnames(gene_range)[1]),
      main = paste("Gene:", gene)
    )
    dev.off()
    message("Plot device closed.")
    if (file.exists(plot_file)) {
      sz <- tryCatch(file.info(plot_file)$size, error = function(e) NA)
      message("Plot saved successfully (bytes=", as.character(sz), ").")
    } else {
      message("Warning: plot file was not written: ", plot_file)
    }
    plots_created <- plots_created + 1L
    
  }, error = function(e) {
    message("Error while generating plot: ", e$message)
  }) 
}
message("Total IGV plots created: ", plots_created)

# Render the HTML report
message("Rendering the HTML report...")
# Ensure required static assets (CSS/JS/images) are available next to the HTML
assets <- c(
  file.path(utilities.path, "report_qc", "tusco_style.css"),
  file.path(utilities.path, "report_qc", "tusco_script.js"),
  file.path(utilities.path, "report_qc", "script.js"),
  file.path(utilities.path, "report_qc", "howToUse.png")
)
invisible(lapply(assets, function(src) {
  if (file.exists(src)) {
    file.copy(src, file.path(output_directory, basename(src)), overwrite = TRUE)
  }
}))
rmarkdown::render(
  input = file.path(utilities.path, "report_qc", "SQANTI3_TUSCO_Report.Rmd"),
  intermediates_dir = output_directory,
  output_dir = output_directory,
  output_file = html.report.file,
  params = params,
  envir = new.env()
)
