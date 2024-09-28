source("analysis/code/generic_analysis.R")

if (FALSE) {
  # How to install
  devtools::install_github("kgori/sigfit")
}

source("analysis/code/common_utils.R")

library(sigfit)

run_sigfit <- function(mut_type) {
  toolname = "sigfit"
  output_home <- file.path("analysis/raw_output", mut_type, toolname, "syn")
  message(toolname," writing to ", output_home)
  
  one_group_of_spectra = function(spectra, 
                         signatures, 
                         more_args) {
    
    
    
    retval = mcmc_samples_fit = sigfit::fit_signatures(
      counts = t(spectra),
      signatures = t(signatures),
      iter = 3000,
      warmup = 1500,
      chains = 4,
      cores = 4,
      seed = global_random_seed
    )
    
    expplus = sigfit::retrieve_pars(retval,  "exposures")
    # activities = sigfit::retrieve_pars(retval, "activities") Vignette only mentins exposures
    
    proportions = t(expplus$mean)
    final_exp = sweep(proportions, 2, colSums(spectra), FUN = "*")
    return(final_exp)
  }

    run_generic_syn(
      dataset_name = mut_type,
      output_home = output_home,
      attribute_function = one_group_of_spectra,
      more_args = NULL
    )
}
