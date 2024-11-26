# Load libraries
library(tidyverse)

# Funtion reads tsv tpm file, get a input tpm_limit threshold, and 
# nonTargetFilter for a threshold.
preProcessTPM <- function(pathToFile, TPM_limit, nonTargetFilter){
  
  # Read TSV file
  matrixTPM <- read_tsv(pathToFile, show_col_types = FALSE) 
  
  # Remove a special column that was in the data as second last column.
  matrixTPM <- matrixTPM |>
    select(- (ncol(matrixTPM) - 1))
  
  # Checking columns are correct.
  print("Columns in dataset:")
  print(colnames(matrixTPM))
  
  # Filtering pipeline 
  FilteredTPM_Matrix <- matrixTPM |>
    mutate(nonTargets = rowSums(across(-c(sample, sample_simple), 
                                  ~ . > TPM_limit)) / (ncol(matrixTPM) - 2)) |>
    filter(nonTargets <= nonTargetFilter) |>
    select(-nonTargets) |> # IF we want to color volcano plot by intervals of this metric it can be included.
    select(sample, sample_simple, everything())
    
  
  return(FilteredTPM_Matrix)
}

# Running pre process function
preProTPM <- preProcessTPM(pathToFile = "data/processed/tcga_counts_subset.tsv.gz", 
                           TPM_limit = 3,
                           nonTargetFilter = 0.6)

# write out data-output from function.
#write_tsv(preProTPM, "data/processed/pre_processed_tpm.tsv.gz")

# collect sample ids needs to subset all our datasets.
sample_to_subset_gtex <- preProTPM |> 
  select(sample)

# read in gtex tpm counts
gtexTPM <- read_tsv("data/raw/filtered_TPM_GTEx.gz")

subsetGtexTPM <- gtexTPM |> 
  semi_join(sample_to_subset_gtex,
            by = c("sample"="sample"))
