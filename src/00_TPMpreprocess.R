library(tidyverse)

### Subsetting using deeploc results
cell_membrane_proteins <- read_tsv(file = ("data/raw/cell_membrane_proteins_enst.txt"))

tcga_counts <- read_tsv(file = ("data/raw/filtered_counts_TCGA.gz"))
enst_list <- cell_membrane_proteins |> select(enst)

tcga_counts_m <- tcga_counts |> 
  mutate(sample_simple = str_remove(sample, "\\.\\d+")) 


tcga_counts_subset <- tcga_counts_m |>
  inner_join(, y = enst_list, by = c("sample_simple" = "enst"))

### Removing rows with column sum below 10
tcga_counts_subset_2 <- tcga_counts_subset |> 
  rowwise() |>
  filter(sum(c_across(where(is.numeric))) >= 10) |>
  ungroup()

### Writing file to disk
tcga_counts_subset_2 |> write_tsv(file = "data/processed/tcga_counts_subset.tsv.gz")


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
write_tsv(preProTPM, "data/processed/pre_processed_tpm.tsv.gz")

# collect sample ids needs to subset all our datasets.
sample_to_subset_gtex <- preProTPM |> 
  select(sample)

# read in gtex tpm counts
gtexTPM <- read_tsv("data/raw/filtered_TPM_GTEx.gz")

subsetGtexTPM <- gtexTPM |> 
  semi_join(sample_to_subset_gtex,
            by = c("sample"="sample"))
