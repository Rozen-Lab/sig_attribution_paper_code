source("analysis/code/generic_analysis.R")

# Only designed to work on SBS
# On data for this project generates error and stops

if (FALSE) {
  # How to install
  devtools::install_github("evenrus/mmsig")
}

source("analysis/code/common_utils.R")

library(mmsig)

run_mmsig <- function(mut_type) {
  toolname = "mmsig"
  output_home <- file.path("analysis/raw_output", mut_type, toolname, "syn")
  message(toolname," writing to ", output_home)
  
  one_group_of_spectra = function(spectra, 
                                  signatures, 
                                  more_args) {
    
    # browser()
    
    spectra = as.data.frame(spectra)
    signatures = as.data.frame(signatures)
    rn= rownames(signatures)
    newcols = data.frame(sub = paste0(substr(rn, 2, 2), ">", substr(rn, 4, 4)),
                         tri = substr(rn, 1, 3))
    signatures = cbind(newcols, signatures)
    
    retval = mmsig::mm_fit_signatures(
      muts.input = spectra,
      input.format = "classes",
      sig.input = signatures,
      sample.sigt.profs = NULL,
      bootstrap = NULL,
      refcheck = FALSE,
      cos_sim_threshold = 0.01,
      force_include = c(),
      dbg = TRUE)
    
    
    browser()
    
    return(exposure_matrix)
  }

    run_generic_syn(
      dataset_name = mut_type,
      output_home = output_home,
      attribute_function = one_group_of_spectra,
      more_args = NULL
    )
}
