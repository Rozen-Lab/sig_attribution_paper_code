# BiocManager::install("YAPSA")

library(YAPSA)

source("analysis/code/generic_analysis.R")

call_yapsa <- function(spectra, 
                       signatures, 
                       more_args
) {
  exposure_matrix <- 
    YAPSA::LCD(in_mutation_catalogue_df = spectra, 
               in_signatures_df = signatures, 
               in_per_sample_cutoff = more_args$in_per_sample_cutoff)
  return(exposure_matrix)
}

run_yapsa <- function(mut_type, in_per_sample_cutoff) {
  # out_fragment = paste0("yapsa_", formatC(in_per_sample_cutoff, digits = 2, format = "f"))
  out_fragment = "yapsa"
  more_args <- list(in_per_sample_cutoff = in_per_sample_cutoff)
  output_home <-
    file.path("analysis/raw_output", mut_type, out_fragment, "syn")
  message("YAPSA, writing to ", output_home)
  
  run_generic_syn(
    dataset_name = mut_type,
    output_home = output_home,
    attribute_function = call_yapsa,
    more_args = more_args
  )
}
