library(tidyverse)
library(here)

# find GTEX TPM limit
plot_gene_expression_density <- function(dataframe, xintercept) {
  dataframe <- dataframe %>%
    mutate(row_mean = 
             (rowMeans(select(., where(is.numeric)))))
  
  dataframe |> ggplot(aes(x = row_mean)) +
    geom_density() +
    xlim(c(-10,100)) + 
    geom_vline(xintercept = xintercept, color = "blue") +
    labs(x = "Mean Gene Expression",
         y = "Density",
         title = "Mean Gene Expression Density")
}

# Change file and xintercept
gtex_tpm <- read_tsv(file = here("data/processed/GTEX_tpm_harmonized.tsv.gz"))
plot_gene_expression_density(dataframe = gtex_tpm, 
                             xintercept = 1)