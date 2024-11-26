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