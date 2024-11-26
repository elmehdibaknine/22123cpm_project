library("tidyverse")

### Subsetting using deeploc results
cell_membrane_proteins <- read_tsv(file = ("data/raw/cell_membrane_proteins_enst.txt"))

tcga_counts <- read_tsv(file = ("data/raw/filtered_TCGA.gz"))
enst_list <- cell_membrane_proteins |> select(enst)

tcga_counts_m <- tcga_counts |> 
  mutate(sample_simple = str_remove(sample, "\\.\\d+")) 


tcga_counts_subset <- tcga_counts_m |>
  inner_join(, y = enst_list, by = c("sample_simple" = "enst"))

### Removing rows with column sum below 10
tcga_counts_subset_2 <- tcga_counts_subset |> 
  rowwise() |>
  filter(sum(c_across(where(is.numeric))) >= 10) |>
  ungroup()

### Writing file to disk
tcga_counts_subset_2 |> write_tsv(file = "data/processed/tcga_counts_subset.tsv.gz")
