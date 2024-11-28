library(reshape2)
library(tidyverse)
library(gplots)
library(RColorBrewer)
library(viridis)

bispecific_correlation_matrix <- function(path_to_TPM_preprocessed){
  
  dds_top_transcripts <- read_rds(path_to_TPM_preprocessed)
  
  matrix_for_correlation <- dds_top_transcripts |> 
    slice_sample(n = 100) |> # remove when running actual analysis and not just samples
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
GTEx_correlation_matrix <- bispecific_correlation_matrix("data/processed/GTEx_tpm_preprocessed.rds")
TCGA_correlation_matrix <- bispecific_correlation_matrix("data/processed/TCGA_tpm_preprocessed.rds")

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