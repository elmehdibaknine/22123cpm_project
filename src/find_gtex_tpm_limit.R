library(tidyverse)
library(here)

# find GTEX TPM limit
plot_gene_expression_density <- function(dataframe, xintercept) {
  dataframe <- dataframe %>%
    mutate(row_mean = 
             (rowMeans(select(., where(is.numeric)))))
  
  dataframe |> ggplot(aes(x = "", y = log10(row_mean + 0.001))) +
    geom_violin() +
    geom_hline(yintercept = log10(0.25 + 0.001)) +
    theme_minimal() +
    scale_fill_viridis_d(option = "plasma") +
    labs(title = "GTEx transcripts mean counts with TPM expression limit Cutoff (0.25)", 
         x = "Transcript means across samples", 
         y = "Log10(transcript_mean + 0.001)")
}

# Change file and xintercept
gtex_tpm <- read_tsv(file = here("data/processed/gtex_tpm_corrected.tsv"))
plot_gene_expression_density(dataframe = gtex_tpm, 
                             xintercept = 1)
