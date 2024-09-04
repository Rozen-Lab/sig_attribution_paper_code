# install.packages(install.packages("SignatureEstimation", repo = NULL, method = "source")
# Package downloaded from https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/index.cgi#signatureestimation
# https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/software/signatureestimation/SignatureEstimation.tar.gz
# Version: 1.0.0 in DESCRIPTION file

library(SignatureEstimation)

source("analysis/code/generic_analysis.R")

call_sigest <- function(spectra, 
                        signatures, 
                        more_args
) {
  allretval <- 
    lapply(colnames(spectra),
           function(sample.id) {
             my.m = spectra[ , sample.id, drop = FALSE]
             retval <- SignatureEstimation::decomposeQP(
               m = my.m,
               P = signatures
             )
             if (any(is.na(retval))) {
               message("NA for ", sample.id)
             }
             return(retval * sum(my.m))
           })
  
  # Need to return an exposure matrix, with samples in columns
  # and signatures in rows

  exposure_matrix <- t(do.call(rbind, allretval))
  colnames(exposure_matrix) = colnames(spectra)
  rownames(exposure_matrix) = colnames(signatures)

  return(exposure_matrix)
}

run_sigest <- function(mut_type) {
  
  output_home <- file.path("analysis/raw_output", mut_type, "sigest/syn")
  time_used <- system.time({
    run_generic_syn(
      dataset_name = mut_type,
      output_home = output_home,
      attribute_function = call_sigest,
      more_args = list()
    )
  })
  
  saveRDS(time_used, file = file.path(output_home, "time_used.Rds"))
}
