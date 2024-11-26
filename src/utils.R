## For retrieving aa sequences from ENSTs
get_aa_seqs_from_ENSTs <- function(transcript_ids, batch_size = 100, seqType = "peptide") {
  # Get human ensembl mart
  message("Creating human mart from biomaRt::useMart with hsapiens_gene_ensembl \n")
  ensembl_hsap <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Initialize an empty tibble to store results
  results <- tibble::tibble()
  
  # Split transcript IDs into batches
  id_batches <- split(transcript_ids, ceiling(seq_along(transcript_ids) / batch_size))
  
  # Loop over each batch
  for (i in seq_along(id_batches)) {
    message("Processing batch", i, "of", length(id_batches), "\n")
    
    # Query the current batch
    batch_result <- biomaRt::getSequence(
      id = id_batches[[i]],
      type = "ensembl_transcript_id",
      seqType = seqType,
      mart = ensembl_hsap
    )
    
    # Convert batch result to a tibble and append to results
    results <- bind_rows(results, as_tibble(batch_result))
  }
  return(results)
}

## For writing fastas
write_fasta_from_named_list <- function(sequences, file_name) {
  # Create a named DNAStringSet or AAStringSet object
  seq_set <- Biostrings::AAStringSet(sequences)
  
  # Write the DNAStringSet object to a FASTA file
  Biostrings::writeXStringSet(seq_set, filepath = file_name)
  
  message("FASTA file written to: ", file_name)
}

### For reading gff3 output from deepTMHMM
tidy_read_gff3 <- function(file_path) {
  # Read the file
  lines <- readLines(file_path)
  
  # Initialize an empty list to store parsed data
  data_list <- list()
  protein_id <- NULL
  length <- NULL
  predicted_tmrs <- NULL
  
  # Loop through each line to parse data
  for (line in lines) {
    
    # Check for comment lines with protein information
    if (grepl("^# sp_", line)) {
      if (grepl("Length:", line)) {
        # Extract protein ID and length
        protein_id <- sub("^# (\\S+).*", "\\1", line)
        length <- as.numeric(sub(".*Length: (\\d+)", "\\1", line))
      }
      if (grepl("Number of predicted TMRs:", line)) {
        # Extract the number of predicted TMRs
        predicted_tmrs <- as.numeric(sub(".*Number of predicted TMRs: (\\d+)", "\\1", line))
      }
      
    } else if (grepl("^sp_", line)) {
      # Parse feature line
      fields <- strsplit(line, "\\s+")[[1]]
      feature_type <- fields[2]
      start <- as.numeric(fields[3])
      end <- as.numeric(fields[4])
      
      # Add parsed information to data list
      data_list <- append(data_list, list(data.frame(
        Protein_ID = protein_id,
        Length = length,
        Predicted_TMRs = predicted_tmrs,
        Feature_Type = feature_type,
        Start = start,
        End = end
      )))
      
    } else if (line == "//") {
      # Reset variables after each protein section
      protein_id <- NULL
      length <- NULL
      predicted_tmrs <- NULL
    }
  }
  
  
  
  # Combine all individual data frames in the list into one
  parsed_data <- do.call(rbind, data_list)
  
  return (parsed_data)
}
