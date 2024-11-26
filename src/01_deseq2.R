# TCGA pancancer
# TCGA TARGET and GTEX
# Script to perform DESeq2 
# 



## Load necessary packages 
library(tidyverse)
library(DESeq2)
library(here)

# load data
tcga_counts <- read_rds(here("data/processed/TCGA_counts_preprocessed.rds"))
gtex_counts <- read_rds(here("data/processed/GTEx_counts_preprocessed.rds"))

#### IF NEEDED ####
# Assume `data` is your matrix with genes as rows and samples as columns
tcga_long <- tcga_counts %>%
  select(-x1214) |> 
  pivot_longer(-transcript, names_to = "sample", values_to = "count")

# Summarize maximum expression for each gene
tcga_gene_summary <- tcga_long %>%
  group_by(transcript) %>%
  summarise(MaxExpression = max(count))


# Visualize maximum gene value across sample distribution
max_plot <- tcga_gene_summary %>%
  ggplot(aes(x = "", y = MaxExpression)) +
  geom_violin(fill = "blue", alpha = 0.5) +
  geom_hline(yintercept = 13) +
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
  filter(MaxExpression > 5)

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
res <- results(dds)

## Convert the results to a data frame and add a column indicating if the results is significant or not
res_df <- res |>
  as_tibble(rownames = "transcript") |> 
  mutate(significant = padj < 0.05)

## Create a volcano plot with ggplot (optional: play around here a bit with colors and themes to improve the readability of your plot)
volcano_plot <- ggplot(res_df, 
       mapping = aes(x = log2FoldChange, 
                     y = -log10(padj), 
                     color = significant)) +
  geom_point() +
  scale_color_manual(values = c("grey", "red")) +
  labs(x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-value", 
       title = "Volcano Plot") + 
  theme_minimal()


ggsave(filename = "volcano_plot.png", plot = volcano_plot)
