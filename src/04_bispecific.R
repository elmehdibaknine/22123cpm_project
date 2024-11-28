library(reshape2)
library(tidyverse)
library(gplots)
library(RColorBrewer)
library(viridis)

bispecific_correlation_matrix <- function(path_to_TPM_preprocessed, transcript_list){
  
  dds_top_transcripts <- read_tsv(path_to_TPM_preprocessed) |>
    filter(transcript %in% transcript_list)

  matrix_for_correlation <- dds_top_transcripts |> 
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

# Load in transcript selection as list
transcript_list <- read_tsv(here("data/processed/ectodomain_containing_transcript.tsv")) |>
  pull(transcript)

# Calling the function on both GTEx and TCGA
GTEx_correlation_matrix <- bispecific_correlation_matrix("data/processed/gtex_tpm_preprocessed.tsv", transcript_list)
TCGA_correlation_matrix <- bispecific_correlation_matrix("data/processed/tcga_tpm_preprocessed.tsv", transcript_list)

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
  mutate(pair = map2_chr(Var1, 
                         Var2, 
                         ~paste(sort(c(.x, .y)), 
                                collapse = "-"))) |>  # Create unique pairs
  distinct(pair, .keep_all = TRUE) |>  # Remove duplicates based on unique pairs
  filter(Freq > 0) |>                  # Keep only positive correlations
  arrange(desc(Freq)) |>               # Sort by correlation value (not absolute)
  slice_head(n = 10) |>                # Select the top 10 correlations
  select(!pair)



#### Plots TCGA vs GTEx
# Create a combined TPM expression dataframe
GTEX_TPM_data <- read_tsv("data/processed/gtex_tpm_preprocessed.tsv")
TCGA_TPM_data <- read_tsv("data/processed/tcga_tpm_preprocessed.tsv")

GTEX_TPM_data_long <- GTEX_TPM_data |>
  pivot_longer(
    -transcript,
    names_to = "sample",
    values_to = "expr_tpm"
  ) |>
  select(-sample) |>
  mutate(data_source = "gtex")

TCGA_TPM_data_long <- TCGA_TPM_data |>
  pivot_longer(
    -transcript,
    names_to = "sample",
    values_to = "expr_tpm"
  ) |>
  select(-sample) |>
  mutate(data_source = "tcga")

combined_tpm <- GTEX_TPM_data_long |>
  rbind(TCGA_TPM_data_long)


# Select a single pair of transcripts from the correlations
pair1 <- correlation_df |>
  dplyr::slice(1) |>
  pull(Var1) |>
  as.character()

pair2 <- correlation_df |>
  dplyr::slice(1) |>
  pull(Var2) |>
  as.character()

pair_data <- combined_tpm |>
  group_by(transcript) %>%
  mutate(pair_group = case_when(
    transcript == pair1 ~ "Pair 1",
    transcript == pair2 ~ "Pair 2",
    TRUE ~ NA_character_
  )) |>
  drop_na()

# Create the plot pair (TCGA vs GTEx expressions)
expr_difference_p <- pair_data |>
  ggplot(aes(x = interaction(transcript, data_source), y = log2(expr_tpm), fill = data_source)) +
    geom_boxplot() +
    facet_wrap(~ pair_group, scales = "free_x") +
    scale_fill_manual(values = c("gtex" = "blue", "tcga" = "red")) +
    theme_minimal() +
    labs(
      title = "TPM expression of CAR-T bispecific set of significantly upregulated transcripts",
      x = "Transcript and Data Source",
      y = "log2(TPM)",
      fill = "Data Source"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 12)
    )

expr_difference_p

TCGA_correlation_matrix |>
  as_tibble(rownames = "transcript") |>
  select(transcript, pair1) |>
  filter(transcript == pair2)

GTEx_correlation_matrix |>
  as_tibble(rownames = "transcript") |>
  select(transcript, pair1) |>
  filter(transcript == pair2)
