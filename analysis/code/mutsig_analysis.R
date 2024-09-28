source("analysis/code/generic_analysis.R")

library(mutSignatures)

run_mutsig <- function(mut_type, cutoff) {
  toolname = "mutsig"
  output_home <- file.path("analysis/raw_output", mut_type, toolname, "syn")
  message(toolname, " writing to ", output_home)
  
  
  one_group_of_spectra <- function(spectra, 
                                   signatures, 
                                   more_args) {
    # browser()
    spec = as.data.frame(spectra)
    spec = mutSignatures::as.mutation.counts(spec)
    sigs = as.data.frame(signatures)
    sigs = mutSignatures::as.mutation.signatures(sigs)
    
    retval <- 
      mutSignatures::resolveMutSignatures(
        mutCountData = spec, 
        signFreqData = sigs)
    retval2 = retval$results$count.result
    
    exposure_matrix = mutSignatures::as.data.frame(retval2, transpose = FALSE)
    
    return(exposure_matrix)
  }
  
  run_generic_syn(
    dataset_name = mut_type,
    output_home = output_home,
    attribute_function = one_group_of_spectra,
    more_args = NULL
  )
}
