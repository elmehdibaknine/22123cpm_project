# Load libraries
library(tidyverse)
library(here)


harmonizeGtex <- function(gtex_tpm_path, gtex_count_path, tpm_harmon_out, count_harmon_out) {
  if (file.exists(tpm_harmon_out) && file.exists(count_harmon_out)) {
    return()
  }
  
  gtex_tpm <- read_tsv(gtex_tpm_path) |>
    mutate(transcript = str_remove(sample, "\\.\\d+"),
           .before = sample) |>
    select(-sample)
  
  gtex_count <- read_tsv(gtex_count_path) |>
    mutate(transcript = str_remove(transcript_id, "\\.\\d+"),
           .before = gene_id) |>
    select(-c(transcript_id, gene_id))
  
  count_transcripts <- gtex_count |>
    pull(transcript)
  
  count_samples <- gtex_count |>
    colnames()
  
  # Write new tpm based on counts data
  gtex_tpm |>
    select(any_of(count_samples)) |>
    filter(transcript %in% count_transcripts) |>
    rename(sample = transcript) |>
    write_tsv(here(tpm_harmon_out))

  
  # Get tpm transcripts and samples
  tpm_transcripts <- gtex_tpm |>
    pull(transcript)
  
  tpm_samples <- gtex_tpm  |>
    colnames()
  
  # Write new counts based on tpms data
  gtex_count |>
    select(any_of(tpm_samples)) |>
    filter(transcript %in% tpm_transcripts) |>
    rename(sample = transcript) |>
    write_tsv(here(count_harmon_out))
  
  return()
}

revertToCounts <- function(log2CountData_path, out_path) {
  if (file.exists(out_path)) {
    return()
  }
  
  log2countdata <- read_tsv(log2CountData_path)
  
  countdata <- log2countdata %>%
    mutate(across(
      where(is.numeric), 
      ~ 2^(.x) - 1
    )) |>
    mutate(across(where(is.numeric), ~ ifelse(. < 0, 0, .)))
  
  countdata |>
    write_tsv(out_path)
}

revertToTPM <- function(log2TPMData_path, out_path) {
  if (file.exists(out_path)) {
    return()
  }
  
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

# Funtion reads tsv tpm file, get a input tpm_limit threshold, and 
# nonTargetFilter for a threshold.
filterHealthyExpressed <- function(GtexTPM, TranscriptSubset, TPM_limit, nonTargetFilter) {
  # Read in GTEx
  gtex_tpm = read_tsv(GtexTPM) |>
    janitor::clean_names() |>
    mutate(transcript = str_remove(sample, "\\.\\d+"),
           .before = sample) |>
    select(-sample)
  
  # Subset GTEx
  gtex_subset <- gtex_tpm |>
    filter(transcript %in% TranscriptSubset)
  
  # Filtering pipeline 
  FilteredTPM_Matrix <- gtex_subset %>%
    mutate(num_expressed = 
             rowMeans(select(., where(is.numeric)) > TPM_limit)
           ) |>
    filter(num_expressed < nonTargetFilter) |>
    select(-num_expressed) |>
    select(transcript, everything())
  
  transcript_list_filtered <- FilteredTPM_Matrix |>
    pull(transcript)
    
  return(transcript_list_filtered)
}

filterDeeplocRowsums <- function(cell_mem_proteins_path, TCGA_counts_path, min_count_sum = 10) {
  ### Subsetting using deeploc results
  cell_membrane_proteins <- read_tsv(file = cell_mem_proteins_path)
  enst_list <- cell_membrane_proteins |> select(enst) |> pull()
  
  
  # Load TCGA
  tcga_counts <- read_tsv(file = TCGA_counts_path) |>
    janitor::clean_names()
  
  tcga_counts <- tcga_counts |>
    mutate(transcript = str_remove(sample, "\\.\\d+"),
           .before = sample) |>
    select(-sample)
  
  tcga_counts_subset_deeploc <- tcga_counts |>
    filter(transcript %in% enst_list)
  
  ### Removing rows with column sum below 10
  tcga_counts_subset_deeploc_rowsums <- tcga_counts_subset_deeploc %>%
    filter(rowSums(select(., where(is.numeric))) >= min_count_sum)
  
  ### Writing file to disk
  transcript_list_filtered <- tcga_counts_subset_deeploc_rowsums |>
    select(transcript) |>
    pull()
  
  return(transcript_list_filtered)
  
}

subsetExpressionData <- function(transcript_subset, expression_data_path, out_path) {
  if (file.exists(out_path)) {
    return()
  }
  expression_data <- read_tsv(expression_data_path) |>
    janitor::clean_names() |>
    mutate(transcript = str_remove(sample, "\\.\\d+"),
           .before = sample) |>
    select(-sample) |>
    select(starts_with("transcript|gtex|tcga"))
  
  expression_data |>
    filter(transcript %in% transcript_subset) |>
    write_rds(out_path)
    
}


# RSEM TPM GTEX log2(x + 0.001)
# RSEM COUNT GTEX log2(x + 1)
# RSEM TPM TCGA log2(x + 0.001)
# RSEM COUNT TCGA log2(x + 1)

# Init values
cell_mem_proteins_path <- here("data/raw/cell_membrane_proteins_enst.txt")
gtex_tpm_raw <- here("data/raw/filtered_tpm_GTEx.gz")
gtex_counts_raw <- here("data/raw/filtered_counts_GTEx.gz")
tcga_tpm_raw <- here("data/raw/filtered_tpm_TCGA.gz")
tcga_counts_raw <- here("data/raw/filtered_counts_TCGA.gz")

gtex_tpm_harmon_out <- here("data/processed/GTEX_tpm_harmonized.tsv.gz")
gtex_count_harmon_out <- here("data/processed/GTEX_counts_harmonized.tsv.gz")

gtex_counts_corrected <- here("data/processed/GTEx_counts_corrected.tsv.gz")
gtex_tpm_corrected <- here("data/processed/GTEx_tpm_corrected.tsv.gz")
tcga_counts_corrected <- here("data/processed/TCGA_counts_corrected.tsv.gz")
tcga_tpm_corrected <- here("data/processed/TCGA_tpm_corrected.tsv.gz")

gtex_counts_processed <- here("data/processed/GTEx_counts_preprocessed.rds")
gtex_tpm_processed <- here("data/processed/GTEx_tpm_preprocessed.rds")
tcga_counts_processed <- here("data/processed/TCGA_counts_preprocessed.rds")
tcga_tpm_processed <- here("data/processed/TCGA_tpm_preprocessed.rds")

# Harmonize shit Gtex datasets
harmonizeGtex(gtex_tpm_raw, 
              gtex_counts_raw,
              gtex_tpm_harmon_out,
              gtex_count_harmon_out)

# Get raw tpm matrices from log2(expr_val + 0.001) transformed state of files currently
revertToTPM(gtex_tpm_harmon_out, gtex_tpm_corrected)
revertToTPM(tcga_tpm_raw, tcga_tpm_corrected)

# Get raw count matrices from log2(expr_val + 1) transformed state of files currently
revertToCounts(gtex_count_harmon_out, gtex_counts_corrected)
revertToCounts(tcga_counts_raw, tcga_counts_corrected)

# Filter tcga based on deeploc subcellular location on membrane and row sums of transcript counts > 10
deeploc_rowsums_transcript <- filterDeeplocRowsums(cell_mem_proteins_path, tcga_counts_corrected)


# Filter gtex based on earlier filter and expression levels in healthy populations
final_transcripts <- filterHealthyExpressed(GtexTPM = gtex_tpm_corrected,
                                            TranscriptSubset = deeploc_rowsums_transcript, 
                                            TPM_limit = 1,
                                            nonTargetFilter = 0.20)

# Write gtex
subsetExpressionData(final_transcripts, gtex_counts_corrected, gtex_counts_processed)
subsetExpressionData(final_transcripts, gtex_tpm_corrected,   gtex_tpm_processed)

# Write TCGA
subsetExpressionData(final_transcripts, tcga_counts_raw, tcga_counts_processed)
subsetExpressionData(final_transcripts, tcga_tpm_corrected, tcga_tpm_processed)



gtex_counts_d <- read_rds(gtex_counts_processed)
gtex_tpm_d <- read_rds(gtex_tpm_processed)


gtex_counts_d |>
  select(-starts_with('gtex'))

gtex_tpm_d |>
  select(-starts_with('gtex'))

tcga_counts_d <- read_rds(tcga_counts_processed)
tcga_tpm_d <- read_rds(tcga_tpm_processed)

tcga_counts_d |>
  select(-starts_with('tcga'))

tcga_tpm_d |>
  select(-starts_with('tcga'))
