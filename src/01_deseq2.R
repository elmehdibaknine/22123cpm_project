# TCGA pancancer
# TCGA TARGET and GTEX
# Script to perform DESeq2 
# 



## Load necessary packages 
library(tidyverse)
library(DESeq2)

# load data
brca_data <- read.csv("data/raw/csv")
gtex_data <- read.csv("data/raw/csv2")

### test ###
load("/home/projects/22123/ex3.Rdata")
brca_long <- COAD_counts

#### IF NEEDED ####
# Assume `data` is your matrix with genes as rows and samples as columns
# brca_long <- brca_data %>%
#   rownames_to_column(var = "Gene") %>%  # Convert row names to a column
#   pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")

#### IF NOT ####
brca_long <- brca_data

# Summarize maximum expression for each gene
brca_gene_summary <- brca_long %>%
  group_by(gene) %>%
  summarise(MaxExpression = max(count))


# Visualize maximum gene value across sample distribution
max_plot <- brca_gene_summary %>%
  ggplot(aes(x = "", y = log10(MaxExpression + 1))) +
  geom_violin(fill = "blue", alpha = 0.5) +
  # geom_jitter(width = 0.2, alpha = 0.5, color = "black") +  # Add jittered points
  theme_minimal() +
  labs(title = "Violin Plot of Max Gene Counts with Points",
       x = "",
       y = "Log10(Max Gene Count + 1)")

ggsave(filename = "max_violin.png", plot = max_plot)

# Visualize brca gene expression distribution per sample
density_plot <- brca_long %>%
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
# Set a threshold, e.g., > 100
filtered_genes <- brca_gene_summary %>%
  filter(MaxExpression > 100)

# Subset the original data to keep only the filtered genes
filtered_brca_long <- brca_long %>%
  filter(gene %in% filtered_genes$gene)

# Check dimensions of the filtered data
dim(brca_long)
dim(filtered_brca_long)




## Convert the counts data to a matrix, where the rownames are genes and the column names are samples
pivot <- function(matrix){
  matrix <- matrix |> 
    pivot_wider(names_from = sample,
                values_from = count) |>
  column_to_rownames("gene")
  return(matrix)
}

filtered_brca_wide <- pivot(filtered_brca_long)
gtex_wide <- pivot(gtex_long)

# inner join cancer and control datasets
counts_matrix <- inner_join(gtex_wide, 
                            filtered_brca_wide, 
                            suffix = c("_gtex", "_brca"))

# Create design data table
disease <- counts_matrix |> 
  select(sample) |> 
  mutate(status = ifelse(grepl("tcga", sample, ignore.case = TRUE),
                         "patient", "control"))



## Create a DESeq data set using overall survival (OS) for the design
dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                              colData = disease,
                              design = ~ status)

## Run DESeq2 on the DESeq data set and look at the results (this make take few minutes, so go and get a coffee)
dds <- DESeq(dds)
res <- results(dds)

## Convert the results to a data frame and add a column indicating if the results is significant or not
res_df <- as_data_frame(res) |> 
  mutate(significant = padj < 0.05)

## Create a volcano plot with ggplot (optional: play around here a bit with colors and themes to improve the readability of your plot)
ggplot(res_df, 
       mapping = aes(x = log2FoldChange, 
                     y = -log10(padj), 
                     color = significant)) +
  geom_point() +
  scale_color_manual(values = c("grey", "red")) +
  labs(x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-value", 
       title = "Volcano Plot") + 
  theme_minimal()