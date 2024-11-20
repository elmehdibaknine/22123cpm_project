library(tidyverse)

preProcessTPM <- function(pathToFile, TPM_limit, nonTargetFilter){
  
  matrixTPM <- read_tsv(pathToFile)
  
  FilteredTPM_Matrix <- matrixTPM |> 
    rowwise() |> 
    mutate(nonTargets = mean(c_across(everything()) > TPM_limit)) |> 
    filter(nonTargets <= nonTargetFilter) |> 
    ungroup()
  
  return(FilteredTPM_Matrix)
}

preProcessTPM(pathToFile = "data/tsv.file", 
              TPM_limit = 3,
              nonTargetFilter = 0.6)