# remotes::install_github("raerose01/deconstructSigs") # Version 1.9.0
library(deconstructSigs)

call_deconstruct <- function(spectra, 
                             signatures, 
                             more_args
                             ) {
  spectra <- t(spectra)
  spectra <- spectra / rowSums(spectra) # This software runs on proporitons
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
               signature.cutoff = more_args$signature_cutoff  # The default is 0.06
             )
           })
  
  # Need to return an exposure matrix, with columns samples
  # and rows signatures
  weights <- lapply(allretval, `[[`, "weights")
  exposure_matrix <- t(do.call(rbind, weights))
  return(exposure_matrix)
}
