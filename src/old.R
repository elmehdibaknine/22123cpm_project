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






########
# deeploc_subset.R
#####

### Subsetting using deeploc results
cell_membrane_proteins <- read_tsv(file = ("data/raw/cell_membrane_proteins_enst.txt"))

tcga_counts <- read_tsv(file = ("data/raw/filtered_TCGA.gz"))
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









######
# 02_select_transcripts.R
#####

library(tidyverse)

# From dds object to data frame
ddsDf <- results(dds,
                 tidy = TRUE)

# significant based on adjusted p-val
ddsDFp1 <- ddsDf |>
  drop_na(padj) |> 
  mutate(significant = padj <= 0.05)

# plot1 - to be modified.
p1 <- ggplot(data = ddsDFp1,
             mapping = aes(x = log2FoldChange,
                           y = -log10(padj),
                           color = significant)) + # color by metric
  geom_point() +
  labs(x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-value", 
       title = "Volcano Plot")

p1

# plot2 - to be modified.
p2 <- ggplot(data = ddsDFp2, 
             mapping = aes(x = log2FoldChange, 
                           y = -log10(padj),
                           color = regulated)) +
  geom_point() +
  labs(x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-value", 
       title = "Volcano Plot")

p2

######################################
# Better Volcano plot ##############
######################################
## Convert the results to a data frame and add a column indicating if the results is significant or not
ddsDf <- results(dds,
                 tidy = TRUE)

# Define thresholds
padj_threshold <- 0.05  # Adjusted p-value threshold
log2fc_threshold <- 1   # Log2 fold change threshold

# Add a new column for color coding
tidy_deseq_results <- ddsDf %>%
  mutate(
    regulation = case_when(
      padj < padj_threshold & log2FoldChange > log2fc_threshold ~ "Up-regulated",
      padj < padj_threshold & log2FoldChange < -log2fc_threshold ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

# Create the volcano plot
ggplot(tidy_deseq_results, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
  geom_point(alpha = 0.8, size = 2) +  
  scale_color_manual(values = c(
    "Up-regulated" = "#F8766D",
    "Down-regulated" = "#619CFF",
    "Not significant" = "gray"
  )) +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted p-value",
    color = "Regulation"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    text = element_text(size = 12)
  )
