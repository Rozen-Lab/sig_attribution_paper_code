# devtools::install_github("gersteinlab/siglasso")

source("analysis/code/generic_analysis.R")

library(siglasso)

call_siglasso <- function(spectra, signatures, more_args) {
  debug = FALSE
  if (debug) {
    if (more_args$cancer_type == "Skin-Melanoma" && (ncol(spectra) > 1)) {
      browser()
      spectra = spectra[ , c("Skin-Melanoma::S.14", "Skin-Melanoma::S.21")]
    }
  }
  
  new_spectra = as.data.frame(spectra)
  new_signatures = as.matrix(as.data.frame(signatures))

  if (more_args$use_prior) {
    message("prior")
    prior = 1 - more_args$sig_props
    prior = prior[colnames(signatures)] # Not all used because of msi/non-msi
    proportion_matrix <- 
      siglasso::siglasso(
        sample_spectrum = new_spectra,
        signature = new_signatures, 
        prior = prior,
        plot = FALSE)
    
  } else {
    message("no prior")
    
    proportion_matrix <- 
      siglasso::siglasso(
        sample_spectrum = new_spectra,
        signature = new_signatures, # Hack to remove SBS96catalog class
        plot = FALSE)
  }

  proportion_assigned = colSums(proportion_matrix)
  remainder = 1 - proportion_assigned
  tmp_proportion_matrix = rbind(remainder, proportion_matrix)
  tmp_exp_matrix = sweep(tmp_proportion_matrix, 2, colSums(new_spectra), FUN = "*")
  exposure_matrix = tmp_exp_matrix[-1, , drop = FALSE]
  colnames(exposure_matrix) = colnames(spectra)
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
