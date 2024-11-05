library(mSigTools)

source("analysis/code/generic_analysis.R")

call_QP_null <- function(spectra, 
                    signatures, 
                    more_args
) {

  
  one_spectrum =function(sample.id) {
    
    # browser()
    my.m = spectra[ , sample.id, drop = FALSE]
    
    retval <- mSigTools::optimize_exposure_QP(
      spectrum = my.m,
      signatures = signatures
    )
    if (any(is.na(retval))) {
      message("NA for ", sample.id)
    }
    return(as.matrix(retval))
  }
           
  allretval <- 
    lapply(colnames(spectra), one_spectrum)
  
  # Need to return an exposure matrix, with samples in columns
  # and signatures in rows
  # browser()
  exposure_matrix <- do.call(cbind, allretval)
  colnames(exposure_matrix) = colnames(spectra)

  return(exposure_matrix)
}

run_QP_null <- function(mut_type) {

  output_home <-
    file.path("analysis/raw_output", mut_type, "QP_null/syn")
  
  
  run_generic_syn(
    dataset_name = mut_type,
    output_home = output_home,
    attribute_function = call_QP_null,
    more_args = NULL
  )
}
