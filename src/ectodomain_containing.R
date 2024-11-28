library(here)
source(here("src/utils.R"))
library(tidyverse)

# Read in the gff3 results from deeptmhmm
deeptmhmm_tidy <- tidy_read_gff3(here("data/processed/TMRs.gff3")) |>
  janitor::clean_names()

# Filter transcripts with ectodomain larger than threshold
ectodomain_length_threshold <- 50
transcript_above_threshold <- deeptmhmm_tidy |>
  group_by(protein_id) |>
  filter(feature_type == "outside") |>
  mutate(ectodomain_length = end - start) |>
  filter(ectodomain_length > ectodomain_length_threshold) |>
  ungroup() |>
  distinct(protein_id) |>
  pull(protein_id)

  
deeptmhmm_tidy |>
  filter(protein_id %in% transcript_above_threshold)
