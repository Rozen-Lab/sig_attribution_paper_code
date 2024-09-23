# devtools::install_github("gersteinlab/siglasso")

source("analysis/code/generic_analysis.R")

library(siglasso)

call_siglasso <- function(spectra, signatures, more_args) {
  new_spectra = as.data.frame(spectra)
  colnames(new_spectra) = colnames(spectra) # The names still get messed up inside siglasso
  new_signatures = as.matrix(as.data.frame(signatures))

  if (more_args$use_prior) {
    message("prior")
    prior = 1 - more_args$sig_props
    prior = prior[colnames(signatures)] # Not all used because of msi/non-msi
    exposure_matrix <- 
      siglasso::siglasso(
        sample_spectrum = new_spectra,
        signature = new_signatures, 
        prior = prior,
        plot = FALSE)
    
  } else {
    message("no prior")
    
    exposure_matrix <- 
      siglasso::siglasso(
        sample_spectrum = new_spectra,
        signature = new_signatures, # Hack to remove SBS96catalog class
        plot = FALSE)
  }
  colnames(exposure_matrix) == colnames(spectra) # Names were modified to make them R identifiers
  return(exposure_matrix)
}

run_siglasso <- function(mut_type, more_args) {
  if (more_args$use_prior) {
    lastbit = "siglasso_wprior"
  } else {
    lastbit = "siglasso"
  }
  output_home <- file.path("analysis/raw_output", mut_type, lastbit, "syn")
  message("siglasso, writing to ", output_home)

    run_generic_syn(
      dataset_name = mut_type,
      output_home = output_home,
      attribute_function = call_siglasso,
      more_args = more_args
    )
}
