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

# Read TUSCO GTF file
message("Reading TUSCO GTF file...")
tusco_gtf <- rtracklayer::import(tusco.file)
tusco_gtf_df <- as.data.frame(tusco_gtf)

# Extract gene-level information
annotation_data <- tusco_gtf_df %>%
  filter(type == "gene") %>%
  select(ensembl = ensembl, refseq = refseq, gene_name = gene_name) %>%
  distinct()

rTUSCO <- nrow(annotation_data)

# Read transcript GTF file
message("Reading transcript GTF file...")
transcript_gtf <- rtracklayer::import(transcript_gtf_file)
transcript_gtf_df <- as.data.frame(transcript_gtf)

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
id_summary <- dplyr::count(classification_data, id_type, sort = TRUE)

message("Summary of id_type counts:")
print(id_summary)

# Identify and display the id_type with the highest count
if (nrow(id_summary) > 0) {
  # Sort by count in descending order
  id_summary_ordered <- id_summary %>% arrange(desc(n))
  
  # If the highest count category is "unknown" and there is at least another category,
  # use the second-highest instead
  if (id_summary_ordered$id_type[1] == "unknown" && nrow(id_summary_ordered) > 1) {
    top_id_type <- id_summary_ordered$id_type[2]
    top_id_n <- id_summary_ordered$n[2]
  } else {
    top_id_type <- id_summary_ordered$id_type[1]
    top_id_n <- id_summary_ordered$n[1]
  }
  
  message(
    "The id_type with the highest count is: ",
    top_id_type,
    " with ",
    top_id_n,
    " entries."
  )
} else {
  message("No id_type classifications were made.")
}

if (top_id_type == "ensembl") {
  # Directly modify the 'gene_id' column
  transcript_gtf_df$gene_id <- sub("\\..*", "", transcript_gtf_df$gene_id)
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
    associated_gene = str_remove(associated_gene, "\\.\\d+$"),         # Remove version number from associated_gene
    associated_transcript = str_remove(associated_transcript, "\\.\\d+$") # Remove version number from associated_transcript
  ) %>%
  distinct(isoform, associated_gene, associated_transcript, .keep_all = TRUE) %>%
  arrange(isoform)

# Metrics and definitions for evaluation against TUSCO
message("Defining TUSCO-related transcript mappings...")
# TUSCO_transcripts: Transcripts mapping to TUSCO genes
TUSCO_transcripts <- classification_data_cleaned %>%
  filter(associated_gene %in% annotation_data[[top_id_type]])

# TUSCO_RM: TUSCO_transcripts that match TUSCO genes as Reference Match (RM)
TUSCO_RM <- TUSCO_transcripts %>%
  filter(associated_gene %in% annotation_data[[top_id_type]] & subcategory == "reference_match")

# Define True Positives (TP): TUSCO transcripts identified as Reference Match (RM)
TP <- TUSCO_RM
TP_TUSCO <- unique(TUSCO_RM$associated_gene)

# Define Partial True Positives (PTP): TUSCO transcripts identified as FSM or ISM but not RM
PTP <- TUSCO_transcripts %>%
  filter(structural_category %in% c("full-splice_match", "incomplete-splice_match") & !associated_gene %in% TP$associated_gene)

# Define False Negatives (FN): TUSCO genes that cannot find any FSM or ISM match
FN <- annotation_data %>%
  filter(!(!!sym(top_id_type) %in% TUSCO_transcripts$associated_gene))

# Define False Positives (FP): Transcripts in NIC, NNC, genic, or fusion categories within TUSCO_transcripts
FP <- TUSCO_transcripts %>%
  filter(structural_category %in% c("novel_in_catalog", "novel_not_in_catalog", "genic", "fusion"))

fsm_ism_count <- TUSCO_transcripts %>%
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
  select(associated_gene, category) %>%
  distinct()

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
    big_category = case_when(
      structural_category == "full-splice_match" & subcategory == "Reference match" ~ "TP", 
      
      # b: PTP
      (structural_category == "full-splice_match" & subcategory %in% c("Alternative 3'end", 
                                                                       "Alternative 5'end",
                                                                       "Alternative 3'5'end")) ~ "PTP",
      structural_category == "incomplete-splice_match" ~ "PTP",
      
      # c: FP
      structural_category %in% c("novel_in_catalog","novel_not_in_catalog","genic_intron",
                                 "genic","antisense","fusion","intergenic") ~ "FP",
      
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
# Generate and save plots for each gene
genes_to_plot <- associated_genes_list$associated_gene
options(ucscChromosomeNames=FALSE)

transcript_gtf_df <- within(transcript_gtf_df, {
  gene_id <- ifelse(
    transcript_id %in% classification_data_cleaned$isoform,
    classification_data_cleaned$associated_gene[match(transcript_id, classification_data_cleaned$isoform)],
    gene_id
  )
})
for (gene in genes_to_plot) {
# Fetch TUSCO gene annotations (exons) for the current gene
  tusco_gene <- tusco_gtf_df %>%
    filter((gene_name == gene | ensembl == gene | refseq == gene ) & type == "exon")
  
  # Fetch transcript exons associated with this gene
  gene_transcripts <- transcript_gtf_df %>%
    filter(gene_id == gene & type == "exon")
  
  message(paste("Creating plot for gene:", gene))
  
  tryCatch({
    message("Starting plot generation for gene: ", gene)
    
    # Merge gene_transcripts with classification_data_cleaned using transcript_id and isoform
    gene_transcripts <- gene_transcripts %>%
      left_join(classification_data_cleaned %>% select(isoform, structural_category), 
                by = c("transcript_id" = "isoform"))
    
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
      select(transcript_id, color) %>%
      distinct()
    color_map <- setNames(transcript_colors$color, transcript_colors$transcript_id)
    mcols(sample_transcripts_gr_unlisted)$fill <- color_map[mcols(sample_transcripts_gr_unlisted)$transcript]
    mcols(sample_transcripts_gr_unlisted)$col <- color_map[mcols(sample_transcripts_gr_unlisted)$transcript]
    
    # Create GRangesList for TUSCO reference transcripts
    message("Creating GRangesList for TUSCO reference transcripts...")
    tusco_transcripts_gr <- makeGRangesListFromDataFrame(
      tusco_gene,
      keep.extra.columns = TRUE,
      seqnames.field = "seqnames",
      start.field = "start",
      end.field = "end",
      strand.field = "strand",
      split.field = "ensembl"
    )
    
    tusco_transcripts_gr_unlisted <- unlist(tusco_transcripts_gr, use.names = TRUE)
    mcols(tusco_transcripts_gr_unlisted)$transcript <- rep(
      names(tusco_transcripts_gr),
      elementNROWS(tusco_transcripts_gr)
    )
    
    # Set fill and col to 'darkgrey' for TUSCO
    mcols(tusco_transcripts_gr_unlisted)$fill <- 'darkgrey'
    mcols(tusco_transcripts_gr_unlisted)$col <- 'darkgrey'
    
    # Define the plotting range
    message("Defining the plotting range...")
  gene_range_start <- min(c(gene_transcripts$start, tusco_gene$start))
  gene_range_end <- max(c(gene_transcripts$end, tusco_gene$end))
    gene_range <- GRanges(
      seqnames = unique(c(gene_transcripts$seqnames, tusco_gene$seqnames)),
      ranges = IRanges(start = gene_range_start, end = gene_range_end),
      strand = unique(c(gene_transcripts$strand, tusco_gene$strand))
    )
    
    genome_axis <- GenomeAxisTrack()
    message("Genome axis created.")
    
    # Generate GeneRegionTrack for sample transcripts
    message("Creating GeneRegionTrack for sample transcripts...")
    sample_tracks <- lapply(unique(mcols(sample_transcripts_gr_unlisted)$transcript), function(transcript) {
      # Subset GRanges for the current transcript
      current_transcript_gr <- sample_transcripts_gr_unlisted[mcols(sample_transcripts_gr_unlisted)$transcript == transcript]
      
      # Fetch the color for the current transcript
      current_color <- unique(mcols(current_transcript_gr)$fill)
      
      # Create a GeneRegionTrack for the current transcript
      GeneRegionTrack(
        current_transcript_gr,
        genome = genome.assembly,
        chromosome = as.character(seqnames(gene_range)[1]),
        name = transcript,
        background.title = current_color,
        fill = current_color,
        col = current_color
      )
    })
    
    # Create GeneRegionTrack for TUSCO reference
    message("Creating GeneRegionTrack for TUSCO reference...")
    tusco_track <- GeneRegionTrack(
      tusco_transcripts_gr_unlisted,
      genome = genome.assembly,
      chromosome = as.character(seqnames(gene_range)[1]),
      name = paste0(gene, " TUSCO Reference"),
      background.title = "darkgrey",  # Set the title background color
      fill = "darkgrey",             # Set fill color
      col = "darkgrey"               # Set border color
    )
    # Combine all transcript tracks into a single list of tracks
    sample_tracks <- unlist(sample_tracks, recursive = FALSE)
    
    # Combine all tracks for plotting
    all_tracks <- c(genome_axis, tusco_track, sample_tracks)
    message("All tracks combined for plotting.")
    
    # Plot and save to file
    plot_file <- file.path(plots_dir, paste0(gene, ".png"))
    message("Saving plot to file: ", plot_file)
    png(plot_file, width = 1000, height = 600)
    plotTracks(
      all_tracks,
      from = gene_range_start,
      to = gene_range_end,
      chromosome = as.character(seqnames(gene_range)[1]),
      main = paste("Gene:", gene)
    )
    dev.off()
    message("Plot saved successfully.")
    
  }, error = function(e) {
    message("Error while generating plot: ", e$message)
  }) 
}

# Render the HTML report
message("Rendering the HTML report...")
rmarkdown::render(
  input = file.path(utilities.path, "report_qc", "SQANTI3_TUSCO_Report.Rmd"),
  intermediates_dir = output_directory,
  output_dir = output_directory,
  output_file = html.report.file,
  params = params,
  envir = new.env()
)
