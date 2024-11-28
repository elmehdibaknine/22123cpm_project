library(reshape2)
library(tidyverse)
library(gplots)
library(RColorBrewer)
library(viridis)

bispecific_correlation_matrix <- function(path_to_TPM_preprocessed, path_to_subset){
  
  dds_top_transcripts <- read_tsv(path_to_TPM_preprocessed)
  most_upregulated_transcripts <- read_tsv(path_to_subset)
  
  matrix_for_correlation <- dds_top_transcripts |> 
    slice_sample(n = 100) |> # remove when running actual analysis and not just samples
    inner_join(most_upregulated_transcripts,
               by = c("transcript" = "transcript")) |> 
    column_to_rownames("transcript") |> 
    t()
  
  # Spearman's correlation matrix
  correlation_matrix <- round(cor(matrix_for_correlation, 
                                  use = "complete.obs",
                                  method = "spearman"), 2)
  # Replacing NA's
  correlation_matrix[is.na(correlation_matrix)] <- 0 
  
  # compressed format
  # melted_correlation_matrix <- melt(correlation_matrix)
  
  return(correlation_matrix)
}

# Calling the function on both GTEx and TCGA
GTEx_correlation_matrix <- bispecific_correlation_matrix("data/processed/gtex_tpm_preprocessed.tsv")
TCGA_correlation_matrix <- bispecific_correlation_matrix("data/processed/tcga_tpm_preprocessed.tsv")

# Correlation diff in BC vs. Healthy
correlation_matrix <- TCGA_correlation_matrix - GTEx_correlation_matrix

# Correlation Matrix melted
melted_correlation_matrix <- melt(correlation_matrix)

# Simple ggplot2 tile plot
ggplot(data = melted_correlation_matrix,
       aes(x = Var1, 
           y = Var2, 
           fill = value)) + 
  geom_tile()

# Generate a color palette using viridis
hmcol <- viridis(100)  # Generate 100 colors from the viridis palette

# Example: Using heatmap.2
heatmap.2(correlation_matrix, 
          trace = 'none', 
          col = hmcol,           # Use viridis colors
          dendrogram = "none",   # No dendrogram
          main = "Heatmap with Viridis Palette")

# Top 10 pairs of CAR's
correlation_df <- as.data.frame(as.table(correlation_matrix)) |> 
  filter(Var1 != Var2) |>             # Remove self-correlations
  mutate(Var1 = pmin(as.character(Var1), 
                     as.character(Var2)),
         Var2 = pmax(as.character(Var1), 
                     as.character(Var2))) |>  # Ensure unique pairs by sorting Var1 and Var2
  distinct(Var1, 
           Var2, 
           .keep_all = TRUE) |>   # Remove duplicates (Var2, Var1)
  filter(Freq > 0) |> 
  arrange(desc(abs(Freq))) |>            # Sort by absolute correlation
  slice_head(n = 10)  

# Plots TCGA vs GTEx
GTEX_TPM_data <- read_tsv("data/processed/gtex_tpm_preprocessed.tsv")
TCGA_TPM_data <- read_tsv("data/processed/tcga_tpm_preprocessed.tsv")

gene_list <- read_tsv("data/processed/transcripts.tsv")

GTEX_TPM_data <- GTEX_TPM_data |> 
  mutate(data_from_GTEx = "GTEx") |> 
  filter(transcript %in% gene_list$transcripts) 

TCGA_TPM_data <- TCGA_TPM_data |> 
  mutate(data_from_TCGA = "TCGA") |> 
  filter(transcript %in% gene_list$transcripts)

TPM_data <- full_join(TCGA_TPM_data, 
                      GTEX_TPM_data,
                      by = "transcript")

  select(transcript, data_from, everything()) |>
  pivot_longer(cols = starts_with("tcga"),
               names_to = "sample_id",
               values_to = "expression")
  
