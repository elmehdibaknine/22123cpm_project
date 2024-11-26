harmonizeGtex <- function(gtex_tpm_path, gtex_count_path, tpm_harmon_out, count_harmon_out) {
  if (file.exists(tpm_harmon_out) && file.exists(count_harmon_out)) {
    return()
  }
  
  gtex_tpm <- read_tsv(gtex_tpm_path) |>
    mutate(transcript = str_remove(sample, "\\.\\d+"),
           .before = sample) |>
    select(-sample)
  
  gtex_count <- read_tsv(gtex_count_path) |>
    mutate(transcript = str_remove(transcript_id, "\\.\\d+"),
           .before = gene_id) |>
    select(-c(transcript_id, gene_id))
  
  count_transcripts <- gtex_count |>
    pull(transcript)
  
  count_samples <- gtex_count |>
    colnames()
  
  # Write new tpm based on counts data
  gtex_tpm |>
    select(any_of(count_samples)) |>
    filter(transcript %in% count_transcripts) |>
    rename(sample = transcript) |>
    write_tsv(here(tpm_harmon_out))
  
  
  # Get tpm transcripts and samples
  tpm_transcripts <- gtex_tpm |>
    pull(transcript)
  
  tpm_samples <- gtex_tpm  |>
    colnames()
  
  # Write new counts based on tpms data
  gtex_count |>
    select(any_of(tpm_samples)) |>
    filter(transcript %in% tpm_transcripts) |>
    rename(sample = transcript) |>
    write_tsv(here(count_harmon_out))
  
  return()
}