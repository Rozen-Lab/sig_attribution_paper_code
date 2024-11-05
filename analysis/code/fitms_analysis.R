# remotes::install_github("Nik-Zainal-Group/signature.tools.lib", ref = "v2.4.5")
# Note, as of 2024 09 03, the DESCRIPTION files still has Version: 2.4.4

library(signature.tools.lib)

source("analysis/code/generic_analysis.R")

global_fitms_rare_sig_thresh = c(0.005, 0.001, 0.05, 0.01, 0.03, 0.06, 0.1, 0.2)

call_fitms <- function(spectra, signatures, more_args) {


  rare_sig_threshold = more_args$rare_sig_threshold
  sig_props <- more_args$sig_props # sig_universe[[cancer_type]]

  # Sanity checking  
  stopifnot(length(setdiff(colnames(signatures), names(sig_props))) == 0)

  common_sig_names <- names(sig_props)[sig_props >= rare_sig_threshold]
  
  ## Some signatures in sig_props may have been removed in non-MSI-high tumors
  common_sig_names <- base::intersect(common_sig_names, colnames(signatures))
  ##                                    
  
  # message("Common ", common_sig_names)

  rare_sig_names <- names(sig_props)[sig_props < rare_sig_threshold]
  
  ##
  rare_sig_names <- base::intersect(rare_sig_names, colnames(signatures))
  ##
  
  # message("rare", rare_sig_names)
  
  if (length(common_sig_names) > 0) {
    common_sig <- signatures[, common_sig_names, drop = FALSE]
  } else {
    stop("No common signatures")
  }
  if (length(rare_sig_names) > 0) {
    rare_sig <- signatures[, rare_sig_names, drop = FALSE]
  } else {
    rare_sig = signatures[, rare_sig_names, drop = FALSE]
  }
  
  retval <-
    signature.tools.lib::FitMS(
      catalogues = spectra,
      commonSignatures = common_sig,
      rareSignatures = rare_sig,
      randomSeed = 2024
    )
  exposures <- t(retval$exposures)
  exposures <- exposures[-which(rownames(exposures) == "unassigned"), , drop = FALSE]
  return(exposures)
}

run_fitms <- function(mut_type, rare_sig_threshold) {
  out_fragment = paste0("fitms_", formatC(rare_sig_threshold, digits = 3, format = "f"))
  more_args <- list(rare_sig_threshold = rare_sig_threshold)
  output_home <-
    file.path("analysis/raw_output", mut_type, out_fragment, "syn")
  message("Writing fitms with rare_sig_threshold ", rare_sig_threshold, " to ", output_home)
  
  run_generic_syn(
    dataset_name = mut_type,
    output_home = output_home,
    attribute_function = call_fitms,
    more_args = more_args
  )
}

