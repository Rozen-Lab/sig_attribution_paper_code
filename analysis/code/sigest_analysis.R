# install.packages(install.packages("SignatureEstimation", repo = NULL, method = "source")
# Package downloaded from https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/index.cgi#signatureestimation
# https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/software/signatureestimation/SignatureEstimation.tar.gz
# Version: 1.0.0 in DESCRIPTION file

library(SignatureEstimation)

source("analysis/code/generic_analysis.R")

call_sigest <- function(spectra, 
                        signatures, 
                        more_args) {

  one_sample = function(sample.id) {
    spectrum = spectra[ , sample.id, drop = FALSE]
    my.m = spectrum / sum(spectrum)

    retval <- SignatureEstimation::decomposeQP(
      m = my.m,
      P = signatures
    )
    if (any(is.na(retval))) {
      message("NA for ", sample.id)
    }
    retval[retval < 1e-8] = 0
    return(retval * sum(spectrum))
  }
  
  allretval <- lapply(colnames(spectra), one_sample)
  
  # Need to return an exposure matrix, with samples in columns
  # and signatures in rows

  exposure_matrix <- t(do.call(rbind, allretval))
  colnames(exposure_matrix) = colnames(spectra)
  rownames(exposure_matrix) = colnames(signatures)

  return(exposure_matrix)
}

run_sigest <- function(mut_type) {
  
  output_home <- file.path("analysis/raw_output", mut_type, "sigest/syn")
  
  run_generic_syn(
    dataset_name = mut_type,
    output_home = output_home,
    attribute_function = call_sigest,
    more_args = list()
  )
}
