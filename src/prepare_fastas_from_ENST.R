library(here)
source(here("src/utils.R"))
library(tidyverse)


ENST_list <- read_tsv(here("data/processed/top_100_transcript.csv")) |>
  pull(transcript)
aa_seqs <- get_aa_seqs_from_ENSTs(ENST_list)

sequences <- aa_seqs |>
  pull(peptide)

names(sequences) <- aa_seqs |>
  pull(ensembl_transcript_id)

write_fasta_from_named_list(sequences, here("data/processed/top_100_transcript.fa"))
