#####################################
#####       BUGSI Report       ######
#####################################
### Author: Tianyuan Liu
### Last Modified: 11/07/2024 by tianyuan.liu@csic.es
#********************** Taking arguments from python script
args <- commandArgs(trailingOnly = TRUE)
class.file <- args[1]
bugsi.file <- args[2]
utilities.path <- args[3]

# Define file paths
# class.file <- '/media/tian/ubuntu/SQANTI_BUGSI/WTC11_cdna_ont_ls/WTC11_cdna_ont_ls_classification.txt'
# bugsi.file <- '/media/tian/ubuntu/GitHub/SQANTI3/utilities/report_qc/bugsi_human.txt'
# utilities.path <- '/media/tian/ubuntu/GitHub/SQANTI3/utilities/'

report.prefix <- strsplit(class.file, "_classification.txt")[[1]][1];
output_directory <- dirname(class.file)
output_name <- basename(report.prefix)
html.report.file <- paste0(output_name, "_BUGSI_report.html");

# Load necessary libraries
message("Loading necessary libraries...")
suppressPackageStartupMessages({
  library(rmarkdown)
  library(tidyr)
  library(dplyr)
  library(readr)
  library(stringr)
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
annotation_data <- read_tsv_safe(
  bugsi.file,
  col_names = c("ensembl", "refseq", "gene_name")
)

# Define regex patterns for ID classification
message("Defining regex patterns for ID classification...")
patterns <- list(
  ensembl = "^(ENSG|ENSMUSG)\\d{11}(\\.\\d+)?$",
  refseq = "^(NM_|NR_|NP_)\\d{6,}$",
  gene_name = "^[A-Z0-9]+$"
)

# Classify associated_gene based on patterns
message("Classifying associated_gene IDs based on their formats...")
classification_data <- classification_data %>%
  mutate(
    id_type = case_when(
      str_detect(associated_gene, patterns$ensembl) ~ "ensembl",
      str_detect(associated_gene, patterns$refseq) ~ "refseq",
      str_detect(associated_gene, patterns$gene_name) ~ "gene_name",
      TRUE ~ "unknown"
    )
  )

# Generate and display summary of id_type counts
id_summary <- classification_data %>%
  count(id_type, sort = TRUE)

message("Summary of id_type counts:")
print(id_summary)

# Identify and display the id_type with the highest count
if (nrow(id_summary) > 0) {
  top_id <- id_summary %>% slice_max(n, n = 1)
  top_id_type <- top_id$id_type
  message(
    "The id_type with the highest count is: ",
    top_id_type,
    " with ",
    top_id$n,
    " entries."
  )
} else {
  message("No id_type classifications were made.")
}

# Clean the classification data
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
  arrange(isoform)

# FSM: full-splice_match
# RM: reference_match
# ISM: incomplete-splice_match
# NIC: novel_in_catalog
# NNC: novel_not_in_catalog
# genic
# fusion

# Metrics and definitions for evaluation against BUGSI
# annotation_data: Ground truth BUGSI model
# BGUSI_transcripts: Transcripts mapping to rBUGSI
# BUGSI_RM: BUGSI_transcripts matching a rBUGSI as Reference Match
# True Positive detections (TP): rBUGSIs identified as RM
# Partial True Positive detections (PTP): rBUGSIs identified as ISM or FSM_non_RM
# False Negative (FN): rBUGSIs can't find any FSM or ISM
# False Positive (FP): NIC + NNC + genic + fusion BUGSI_transcripts

# Non_redundant Precision: TP/BUGSI_transcripts
# Positive Detection Rate: unique(TP+PTP)/rBUGSIs
# False Discovery Rate: (BUGSI_transcripts - BUGSI_RM)/BUGSI_transcripts
# False Detection Rate: FP/BUGSI_transcripts
# Redundancy: (FSM + ISM)/unique(TP+PTP)

message("Defining BUGSI-related transcript mappings...")
# BUGSI_transcripts: Transcripts mapping to BUGSI genes
BUGSI_transcripts <- classification_data_cleaned %>%
  filter(associated_gene %in% annotation_data[[top_id_type]])

# BUGSI_RM: BUGSI_transcripts that match BUGSI genes as Reference Match (RM)
BUGSI_RM <- BUGSI_transcripts %>%
  filter(associated_gene %in% annotation_data[[top_id_type]] & subcategory == "reference_match")

# Define True Positives (TP): BUGSI transcripts identified as Reference Match (RM)
TP <- BUGSI_RM

# Define Partial True Positives (PTP): BUGSI transcripts identified as FSM or ISM but not RM
PTP <- BUGSI_transcripts %>%
  filter(structural_category %in% c("full-splice_match", "incomplete-splice_match") & !associated_gene %in% TP$associated_gene)

# Define False Negatives (FN): BUGSI genes that cannot find any FSM or ISM match
FN <- annotation_data %>%
  filter(!(!!sym(top_id_type) %in% BUGSI_transcripts$associated_gene))

# Define False Positives (FP): Transcripts in NIC, NNC, genic, or fusion categories within BUGSI_transcripts
FP <- BUGSI_transcripts %>%
  filter(structural_category %in% c("novel_in_catalog", "novel_not_in_catalog", "genic", "fusion"))

fsm_ism_count <- BUGSI_transcripts %>%
  filter(structural_category %in% c("full-splice_match", "incomplete-splice_match")) %>%
  nrow()

# Remove the specified columns and any columns that are all NA
TP <- TP %>%
  select(-associated_transcript, -id_type) %>%
  select_if(~ !all(is.na(.)))

PTP <- PTP %>%
  select(-associated_transcript, -id_type) %>%
  select_if(~ !all(is.na(.)))

FP <- FP %>%
  select(-associated_transcript, -id_type) %>%
  select_if(~ !all(is.na(.)))

# Calculate metrics
message("Calculating evaluation metrics...")


# Initialize variables to store metric values
non_redundant_precision <- NA
positive_detection_rate <- NA
false_discovery_rate <- NA
false_detection_rate <- NA
redundancy <- NA

# Calculate Non-redundant Precision: TP / BUGSI_transcripts
if (nrow(BUGSI_transcripts) > 0) {
  non_redundant_precision <- nrow(TP) / nrow(BUGSI_transcripts)
} else {
  warning("BUGSI_transcripts has zero rows. Non-redundant Precision set to NA.")
}

# Calculate Precision: (TP + PTP) / BUGSI_transcripts
message("Calculating Precision...")
if (fsm_ism_count > 0) {
  redundant_precision <- (nrow(TP) + nrow(PTP)) /  nrow(BUGSI_transcripts)
  message("Precision calculated successfully.")
} else {
  redundant_precision <- NA
  warning("FSM + ISM count is zero. Redundant Precision set to NA.")
}

# Calculate Positive Detection Rate: unique(TP + PTP) / rBUGSIs
total_rBUGSIs <- nrow(annotation_data)
unique_detected_genes <- length(unique(c(TP$associated_gene, PTP$associated_gene)))
if (total_rBUGSIs > 0) {
  positive_detection_rate <- unique_detected_genes / total_rBUGSIs
} else {
  warning("annotation_data has zero rows. Positive Detection Rate set to NA.")
}

# Calculate False Discovery Rate: (BUGSI_transcripts - BUGSI_RM) / BUGSI_transcripts
if (nrow(BUGSI_transcripts) > 0) {
  false_discovery_rate <- (nrow(BUGSI_transcripts) - nrow(BUGSI_RM)) / nrow(BUGSI_transcripts)
} else {
  warning("BUGSI_transcripts has zero rows. False Discovery Rate set to NA.")
}

# Calculate False Detection Rate: FP / BUGSI_transcripts
if (nrow(BUGSI_transcripts) > 0) {
  false_detection_rate <- nrow(FP) / nrow(BUGSI_transcripts)
} else {
  warning("BUGSI_transcripts has zero rows. False Detection Rate set to NA.")
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
    "Non-redundant Precision",
    "Precision",
    "Positive Detection Rate",
    "False Discovery Rate",
    "False Detection Rate",
    "Redundancy"
  ),
  Value = c(
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

# Render the HTML report
message("Rendering the HTML report...")
rmarkdown::render(
  input = paste(utilities.path, "/report_qc/SQANTI3_BUGSI_Report.Rmd", sep = "/"),
  intermediates_dir = output_directory,
  output_dir = output_directory,
  output_file = html.report.file
)
