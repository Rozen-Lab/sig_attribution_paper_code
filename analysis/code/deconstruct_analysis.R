# remotes::install_github("raerose01/deconstructSigs") # Version 1.9.0
library(deconstructSigs)

source("analysis/code/generic_analysis.R")

call_deconstruct <- function(spectra, 
                             signatures, 
                             more_args
                             ) {
  spectra <- t(spectra)
  spectra <- spectra / rowSums(spectra) # This software runs on proportions
  spectra <- as.data.frame(spectra)
  signatures <-as.data.frame(t(signatures))

  allretval <- 
    lapply(rownames(spectra),
           function(sample.id) {
             retval <- deconstructSigs::whichSignatures(
               tumor.ref = spectra,
               sample.id = sample.id,
               signatures.ref = signatures,
               contexts.needed = FALSE, # Not counts, proportions in the spectra
               tri.counts.method = "default", # No adjustments for exome versus genome 
               # trinucleotide frequencies
               signature.cutoff = more_args$signature.cutoff  # The default is 0.06
             )
           })
  
  # Need to return an exposure matrix, with columns samples
  # and rows signatures
  weights <- lapply(allretval, `[[`, "weights")
  exposure_matrix <- t(do.call(rbind, weights))
  exposure_matrix <- exposure_matrix * sum(spectra)
  return(exposure_matrix)
}

run_deconstruct <- function(mut_type, signature.cutoff) {
  out_fragment = paste0("deconstruct_", 
                        formatC(signature.cutoff, digits = 2, format = "f"))
  more_args = list(signature.cutoff = signature.cutoff) 
  output_home <-
    file.path("analysis/raw_output", mut_type, out_fragment, "syn")
  message("deconstructSigs, writing to ", output_home)
  
  run_generic_syn(
    dataset_name = mut_type,
    output_home = output_home,
    attribute_function = call_deconstruct,
    more_args = more_args
  )
}
