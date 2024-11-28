# TCGA pancancer
# TCGA TARGET and GTEX
# Script to perform DESeq2 
# 



## Load necessary packages 
library(tidyverse)
library(DESeq2)
library(here)

# Load counts data
tcga_counts <- read_tsv(here("data/processed/tcga_counts_preprocessed.tsv"))
gtex_counts <- read_tsv(here("data/processed/gtex_counts_preprocessed.tsv"))


#### IF NEEDED ####
# Assume `data` is your matrix with genes as rows and samples as columns
tcga_long <- tcga_counts %>%
  pivot_longer(-transcript, names_to = "sample", values_to = "count")

# Summarize maximum expression for each gene
tcga_gene_summary <- tcga_long %>%
  group_by(transcript) %>%
  summarise(MaxExpression = max(count))


# Visualize maximum gene value across sample distribution
max_plot <- tcga_gene_summary %>%
  ggplot(aes(x = "", y = log10(MaxExpression +1))) +
  geom_violin(fill = "blue", alpha = 0.5) +
  geom_hline(yintercept = log10(10)) +
  # geom_jitter(width = 0.2, alpha = 0.5, color = "black") +  # Add jittered points
  theme_minimal() +
  labs(title = "Violin Plot of Max Gene Counts with Points",
       x = "",
       y = "Max Gene Count")

ggsave(filename = "max_violin.png", plot = max_plot)

# Visualize tcga gene expression distribution per sample
density_plot <- tcga_long %>%
  ggplot(aes(x = log10(count + 1), color = sample)) +
  geom_density(alpha = 0.5) +
  # scale_x_log10() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Density Distribution of Gene Counts",
       x = "Gene Counts",
       y = "Density",
       color = "Sample")

ggsave(filename = "density_plot.png", plot = density_plot)


## filter genes that are expressed lower than a threshhold across all samples 
# Set a threshold, e.g., > 5
filtered_genes <- tcga_gene_summary %>%
  filter(MaxExpression > 10)

# Subset the original data to keep only the filtered genes
filtered_tcga_long <- tcga_long %>%
  filter(transcript %in% filtered_genes$transcript)

# Check dimensions of the filtered data
dim(tcga_long)
dim(filtered_tcga_long)




## Convert the counts data to a matrix, where the rownames are genes and the column names are samples
pivot <- function(matrix){
  matrix <- matrix |> 
    pivot_wider(names_from = sample,
                values_from = count)
  return(matrix)
}

filtered_tcga_wide <- pivot(filtered_tcga_long)
gtex_wide <- gtex_counts

# inner join cancer and control datasets
counts_matrix <- inner_join(gtex_wide, 
                            filtered_tcga_wide,
                            by = "transcript")

# Create design data table
disease <- counts_matrix |> 
  select(-transcript) |> 
  colnames() |>  
  as.data.frame() 


disease <- disease |> 
  mutate(sample = colnames(select(counts_matrix, -transcript)),
         status = ifelse(grepl("tcga", sample, ignore.case = TRUE),
                         "patient", "control")) |> 
  select(sample, status)

# 
counts_matrix <- counts_matrix |> 
  column_to_rownames("transcript") 

# round counts to integers
counts_matrix <- ceiling(counts_matrix)
## Create a DESeq data set using health status (status) for the design
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = disease,
                              design = ~ status)

## Run DESeq2 on the DESeq data set and look at the results (this make take few minutes, so go and get a coffee)
dds <- DESeq(dds)
## Convert to tidy dataframe
res_tidy <- results(dds, tidy = TRUE)


# Define thresholds
padj_threshold <- 0.05  # Adjusted p-value threshold
log2fc_threshold <- 1   # Log2 fold change threshold

# Add a new column for color coding
tidy_deseq_results <- res_tidy %>%
  mutate(
    regulation = case_when(
      padj < padj_threshold & log2FoldChange > log2fc_threshold ~ "Up-regulated",
      padj < padj_threshold & log2FoldChange < -log2fc_threshold ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

# Create the volcano plot
volcano <- ggplot(tidy_deseq_results, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
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


ggsave(filename = "volcano.png", plot = volcano)


# Extract top 100
top_transcripts <- res_tidy |> 
  dplyr::rename(transcript = row) |>
  mutate(metric = -log10(padj) * log2FoldChange) |> 
  arrange(desc(metric)) |>
  head(n = 100) |>
  select(transcript)

top_transcripts |>
  write_tsv(file = here("data/processed/top_upregulated_transcripts.tsv"))
