library(here)
source(here("src/utils.R"))
library(tidyverse)

# Read in the gff3 results from deeptmhmm
deeptmhmm_tidy <- tidy_read_gff3(here("data/processed/TMRs.gff3")) |>
  janitor::clean_names() |>
  dplyr::rename(transcript = protein_id)

# Filter transcripts with ectodomain larger than threshold
ectodomain_length_threshold <- 50
deeptmhmm_tidy_above_threshold <- deeptmhmm_tidy |>
  group_by(transcript) |>
  filter(feature_type == "outside") |>
  mutate(ectodomain_length = end - start) |>
  filter(ectodomain_length > ectodomain_length_threshold) |>
  ungroup() |>
  distinct(transcript)

deeptmhmm_tidy_above_threshold |>
  write_tsv(here("data/processed/ectodomain_containing_transcript.tsv"))

# Pull a list of protein ids above threshold
transcript_above_threshold <- deeptmhmm_tidy_above_threshold |>
  pull(transcript)

# Plot selected transcripts predicted topology features
deeptmhmm_tidy |>
  filter(transcript %in% transcript_above_threshold[1:9]) |>
  ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = feature_type)) +
  geom_rect() +
  facet_wrap(~ transcript, scales = "free_y") +   # Facet by transcript (transcript)
  scale_fill_manual(values = c("signal" = "orange", "outside" = "green", "TMhelix" = "darkgrey", "inside" = "red")) + # Customize colors
  theme_minimal() +
  theme(
    axis.text.y = element_blank(), # Hide y-axis text (only horizontal bars)
    axis.ticks.y = element_blank(), # Hide y-axis ticks
    strip.text.x = element_text(angle = 0, hjust = 0.5), # Rotate facet labels if needed
    panel.grid = element_blank() # Remove grid lines
  ) +
  labs(x = "Position", y = NULL, fill = "Feature Type") +
  ggtitle("Protein Feature Types by Transcript")
