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


## Convert the counts data to a matrix, where the rownames are genes and the column names are samples
counts_matrix <- counts_matrix |> 
  pivot_wider(names_from = sample,
              values_from = count) |> 
  column_to_rownames("gene")

## Create a DESeq data set using overall survival (OS) for the design
dds <- DESeqDataSetFromMatrix(countData = COAD_counts_matrix,
                              colData = surv_filtered,
                              design = ~ OS)

## Run DESeq2 on the DESeq data set and look at the results (this make take few minutes, so go and get a coffee)
dds <- DESeq(dds)
res <- results(dds)

## Convert the results to a data frame and add a column indicating if the results is significant or not
res_df <- as_data_frame(res) |> 
  mutate(significant = padj < 0.05)