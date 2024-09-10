library(mSigTools)

source("analysis/code/generic_analysis.R")

call_QP <- function(spectra, 
                    signatures, 
                    more_args
) {

  allretval <- 
    lapply(colnames(spectra),
           function(sample.id) {
             
             # browser()
             my.m = spectra[ , sample.id, drop = FALSE]
             
             if (more_args$trim_less_than_is_fraction) {
               more_args$trim_less_than = sum(my.m) * more_args$trim_less_than
             }
             
             retval <- mSigTools::find_best_reconstruction_QP(
               target.sig = my.m,
               sig.universe = signatures,
               trim.less.than = more_args$trim_less_than
             )
             if (any(is.na(retval))) {
               message("NA for ", sample.id)
             }
             return(retval$optimized.exposure)
           })
  
  # Need to return an exposure matrix, with samples in columns
  # and signatures in rows
  # browser()
  exposure_matrix <- t(do.call(rbind, allretval))
  colnames(exposure_matrix) = colnames(spectra)

  return(exposure_matrix)
}

run_QP <- function(mut_type, QP_more_args) {
  
  out_fragment = 
    paste0("QP_", 
           formatC(QP_more_args$trim_less_than, digits = 2, format = "f"))
  output_home <-
    file.path("analysis/raw_output", mut_type, out_fragment, "syn")
  
  
  run_generic_syn(
    dataset_name = mut_type,
    output_home = output_home,
    attribute_function = call_QP,
    more_args = QP_more_args
  )
}
