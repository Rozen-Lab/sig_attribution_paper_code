source("analysis/code/generic_analysis.R")

source("analysis/code/common_utils.R")

library(SigsPack)

run_sigspack <- function(mut_type) {
  toolname = "sigspack"
  output_home <- file.path("analysis/raw_output", mut_type, toolname, "syn")
  message(toolname," writing to ", output_home)
  
  one_group_of_spectra = function(spectra, 
                                  signatures, 
                                  more_args) {
    # browser()
    retval = SigsPack::signature_exposure(
      mut_cat = spectra,
      P = signatures)
    # browser()
    exposures = retval$exposures
    colnames(exposures) = colnames(spectra)
    final_exp = sweep(exposures, 2, colSums(spectra), FUN = "*")
    return(final_exp)
  }
  
  run_generic_syn(
    dataset_name = mut_type,
      output_home = output_home,
      attribute_function = one_group_of_spectra,
      more_args = NULL
    )
}
