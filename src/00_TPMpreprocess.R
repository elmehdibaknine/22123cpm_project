# Load libraries
library(tidyverse)
library(here)

options(readr.show_progress = FALSE, readr.show_col_types = FALSE)

revertToCounts <- function(log2CountData_path, out_path) {
  
  log2countdata <- read_tsv(log2CountData_path)
  
  countdata <- log2countdata %>%
    mutate(across(
      where(is.numeric), 
      ~ 2^(.x) - 1
    )) |>
    mutate(across(where(is.numeric), ~ ifelse(. < 0, 0, .))) |>
    mutate(across(where(is.numeric), round))
  
  countdata |>
    write_tsv(out_path)
}

revertToTPM <- function(log2TPMData_path, out_path) {
  
  log2TPMdata <- read_tsv(log2TPMData_path)
  
  countdata <- log2TPMdata %>%
    mutate(across(
      where(is.numeric), 
      ~ 2^(.x) - 0.001
    )) |>
    mutate(across(where(is.numeric), ~ ifelse(. < 0, 0, .)))
  
  countdata |>
    write_tsv(out_path)
}

filterRowsums <- function(TCGA_counts_path, GTEX_counts_path, min_count_sum = 10) {
  # Load Counts
  tcga_counts <- read_tsv(file = TCGA_counts_path) |>
    janitor::clean_names() |>
    rename(transcript = sample)
  
  gtex_counts <- read_tsv(file = GTEX_counts_path) |>
    janitor::clean_names() |>
    rename(transcript = sample)
  
  ### Removing rows with column sum below threshold
  tcga_counts_subset_rowsums <- tcga_counts %>%
    filter(rowSums(select(., where(is.numeric))) >= min_count_sum)
  
  gtex_counts_subset_rowsums <- gtex_counts %>%
    filter(rowSums(select(., where(is.numeric))) >= min_count_sum)
  
  ### Get transcripts from tcga and gtex
  tcga_transcript_filtered <- tcga_counts_subset_rowsums |>
    pull(transcript)
  
  gtex_transcript_filtered <- gtex_counts_subset_rowsums |>
    pull(transcript)
  
  in_common_transcript_filtered <- intersect(tcga_transcript_filtered, gtex_transcript_filtered)
  
  return(in_common_transcript_filtered)
  
}

# Funtion reads tsv tpm file, get a input tpm limit threshold, and 
filterHealthyExpressed <- function(GtexTPM, TranscriptSubset, TPMLimit, nonTargetFilter) {
  # Read in GTEx
  gtex_tpm = read_tsv(GtexTPM) |>
    janitor::clean_names() |>
    rename(transcript = sample)
  
  # Subset GTEx
  gtex_subset <- gtex_tpm |>
    filter(transcript %in% TranscriptSubset)
  
  # Filtering pipeline 
  FilteredTPM_Matrix <- gtex_subset %>%
    mutate(num_expressed = 
             rowMeans(select(., where(is.numeric)) > TPMLimit)
           ) |>
    filter(num_expressed < nonTargetFilter) |>
    select(-num_expressed) |>
    select(transcript, everything())
  
  transcript_list_filtered <- FilteredTPM_Matrix |>
    pull(transcript)
    
  return(transcript_list_filtered)
}

subsetExpressionData <- function(transcript_subset, expression_data_path, out_path) {
  expression_data <- read_tsv(expression_data_path) |>
    janitor::clean_names() |>
    rename(transcript = sample)
  
  expression_data |>
    filter(transcript %in% transcript_subset) |>
    write_tsv(out_path)
    
}


# RSEM TPM GTEX log2(x + 0.001)
# RSEM COUNT GTEX log2(x + 1)
# RSEM TPM TCGA log2(x + 0.001)
# RSEM COUNT TCGA log2(x + 1)

# Init values
gtex_tpm_raw <-           here("data/processed/gtex_tpm_subset.tsv")
gtex_counts_raw <-        here("data/processed/gtex_counts_subset.tsv")
gtex_counts_corrected <-  here("data/processed/gtex_counts_corrected.tsv")
gtex_tpm_corrected <-     here("data/processed/gtex_tpm_corrected.tsv")
gtex_counts_processed <-  here("data/processed/gtex_counts_preprocessed.tsv")
gtex_tpm_processed <-     here("data/processed/gtex_tpm_preprocessed.tsv")


tcga_tpm_raw <-           here("data/processed/tgca_tpm_subset.tsv")
tcga_counts_raw <-        here("data/processed/tcga_counts_subset.tsv")
tcga_counts_corrected <-  here("data/processed/tcga_counts_corrected.tsv")
tcga_tpm_corrected <-     here("data/processed/tcga_tpm_corrected.tsv")
tcga_counts_processed <-  here("data/processed/tcga_counts_preprocessed.tsv")
tcga_tpm_processed <-     here("data/processed/tcga_tpm_preprocessed.tsv")

# Get raw tpm matrices from log2(expr_val + 0.001) transformed state of files currently
revertToTPM(gtex_tpm_raw, gtex_tpm_corrected)
revertToTPM(tcga_tpm_raw, tcga_tpm_corrected)

# Get raw count matrices from log2(expr_val + 1) transformed state of files currently
revertToCounts(gtex_counts_raw, gtex_counts_corrected)
revertToCounts(tcga_counts_raw, tcga_counts_corrected)

# Filter tcga based on row sums of transcript counts > 10
rowsums_transcript <- filterRowsums(tcga_counts_corrected, gtex_counts_corrected)


# Filter gtex based on earlier filter and expression levels in healthy populations
final_transcripts <- filterHealthyExpressed(GtexTPM = gtex_tpm_corrected,
                                            TranscriptSubset = rowsums_transcript, 
                                            TPMLimit = 0.25,
                                            nonTargetFilter = 0.20)

# Write gtex
subsetExpressionData(final_transcripts, gtex_counts_corrected, gtex_counts_processed)
subsetExpressionData(final_transcripts, gtex_tpm_corrected,   gtex_tpm_processed)

# Write TCGA
subsetExpressionData(final_transcripts, tcga_counts_corrected, tcga_counts_processed)
subsetExpressionData(final_transcripts, tcga_tpm_corrected, tcga_tpm_processed)

gtex_c <- read_tsv(gtex_counts_processed)
