source("analysis/code/generic_analysis.R")

library(mutSignatures)

call_mutsig <- function(spectra, 
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

run_mutsig <- function(mut_type, cutoff) {
  output_home <- file.path("analysis/raw_output", mut_type, "mutsig/syn")
  message("mutSignatures, writing to ", output_home)

    run_generic_syn(
      dataset_name = mut_type,
      output_home = output_home,
      attribute_function = call_mutsig,
      more_args = NULL
    )
}
