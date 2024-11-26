library(here)
source(here("src/utils.R"))
library(tidyverse)


load(here("../weekly_exercises/data/ex9/HER2_isoform_expression.RData"))


peptide_tib <- enst_expr |>
  pull(rowname) |>
  get_aa_seqs_from_ENSTs()

sequences <- peptide_tib |>
  pull(peptide)

names(sequences) <- peptide_tib |>
  pull(ensembl_transcript_id)

write_fasta_from_named_list(sequences, here("data/processed/test.fa"))


tidy_read_gff3(here("data/processed/test.gff3"))
