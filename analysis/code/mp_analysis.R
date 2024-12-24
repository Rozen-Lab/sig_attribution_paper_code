library(MutationalPatterns)

source("analysis/code/generic_analysis.R")

call_mp <- function(spectra, 
                    signatures, 
                    more_args) {
  retval <- 
    MutationalPatterns::fit_to_signatures_strict(
      mut_matrix = spectra, 
      signatures = signatures)
  exposure_matrix <- retval$fit_res$contribution
  return(exposure_matrix)
}

run_mp <- function(mut_type, cutoff) {
  output_home <- file.path("analysis/raw_output", mut_type, "mp/syn")
  message("MutationalPatterns, writing to ", output_home)

    run_generic_syn(
      dataset_name = mut_type,
      output_home = output_home,
      attribute_function = call_mp,
      more_args = NULL
    )
}
