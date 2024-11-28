library(here)
source(here("src/utils.R"))
library(tidyverse)

# Read in the gff3 results from deeptmhmm
tidy_read_gff3(here("data/processed/test.gff3"))