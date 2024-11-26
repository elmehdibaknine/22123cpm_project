# Load libraries
library(tidyverse)
library(here)


harmonizeGtex <- function(gtex_tpm_path, gtex_count_path) {
  tpm_harmon_out <- "data/processed/GTEX_tpm_harmonized.tsv.gz"
  count_harmon_out <- "data/processed/GTEX_counts_harmonized.tsv.gz"
  if (file.exists(here(tpm_harmon_out)) && file.exists(here(count_harmon_out))) {
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
    write_tsv(here(counts_harmon_out))
  
  return()
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

filterDeeplocRowsums <- function(min_count_sum = 10) {
  ### Subsetting using deeploc results
  cell_membrane_proteins <- read_tsv(file = (here("data/raw/cell_membrane_proteins_enst.txt")))
  enst_list <- cell_membrane_proteins |> select(enst) |> pull()
  
  
  # Load TCGA
  tcga_counts <- read_tsv(file = here("data/raw/filtered_counts_TCGA.gz")) |>
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
  expression_data <- read_tsv(expression_data_path) |>
    janitor::clean_names() |>
    mutate(transcript = str_remove(sample, "\\.\\d+"),
           .before = sample) |>
    select(-sample)
  
  expression_data |>
    filter(transcript %in% transcript_subset) |>
    write_rds(out_path)
    
}

# Harmonize shit Gtex datasets
harmonizeGtex(here("data/raw/filtered_tpm_GTEx.gz"), here("data/raw/filtered_counts_GTEx.gz"))

# Filter tcga based on deeploc subcellular location on membrane and row sums of transcript counts > 10
deeploc_rowsums_transcript <- filterDeeplocRowsums()


# Filter gtex based on earlier filter and expression levels in healthy populations
final_transcripts <- filterHealthyExpressed(GtexTPM = "data/processed/GTEx_tpm_harmonized.tsv.gz",
                                            TranscriptSubset = deeploc_rowsums_transcript, 
                                            TPM_limit = 1,
                                            nonTargetFilter = 0.20)

subsetExpressionData(final_transcripts, here("data/processed/GTEx_tpm_harmonized.tsv.gz"),    here("data/processed/GTEx_tpm_preprocessed.rds"))
subsetExpressionData(final_transcripts, here("data/processed/GTEx_counts_harmonized.tsv.gz"), here("data/processed/GTEx_counts_preprocessed.rds"))
subsetExpressionData(final_transcripts, here("data/raw/filtered_tpm_TCGA.gz"),                here("data/processed/TCGA_tpm_preprocessed.rds"))
subsetExpressionData(final_transcripts, here("data/raw/filtered_counts_TCGA.gz"),             here("data/processed/TCGA_counts_preprocessed.rds"))

