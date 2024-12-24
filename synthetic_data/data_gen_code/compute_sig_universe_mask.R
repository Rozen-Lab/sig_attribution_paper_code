# This script computes an "sig_universe_mask.csv" file for each of
# SBS, DBS, and ID. The sampe of sig_universe_mask is the same as
# exposures or all tumor for one too, i.e. <number of signatures> X 900
# The mask contains a 1 in a cell that signature was in the
# signature universe for the sample in that column. We need this because,
# for a given cancer type, MSI samples have more signaures in the universe.

# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())

source("analysis/code/generic_analysis.R")
library(mSigTools)
library(mSigAct)

compute_sig_universe_mask <- function(dataset_name, output_home) {

  synthetic_data_folder = "synthetic_data"
  all_inputs <- get_all_input( # defined in generic_analysis.R
    dataset_name = dataset_name,
    data_top_folder_name = synthetic_data_folder
  )
  all_spectra <- all_inputs$spectra_list
  # gt_sig <- all_inputs$ground_truth_sigs
  sig_universe <- all_inputs$signature_universes
  
  cancer_types <- all_inputs$cancer_types
  
  all_sig_universes = list()

  for (cancer_type in cancer_types) { 
    message(cancer_type)
    sample_names <- colnames(all_spectra[[cancer_type]])
    sig_names <- names(sig_universe[[cancer_type]]) # These are signatures seen in this cancer type

    sig_universes <-
      sort_non_msi_and_msi(
        sample_names = sample_names,
        sig_names = sig_names)

    all_sig_universes = c(all_sig_universes, list(sig_universes))
    
  } # for (cancer_type in cancer_types...
  # browser()
  final_sig_universe = mSigAct:::MergeListOfExposures(all_sig_universes)
  
  original_spectra = all_inputs$all_spectra
  
  if (setequal(colnames(original_spectra), colnames(final_sig_universe))) {
    final_sig_universe <-
      final_sig_universe[, colnames(original_spectra), drop = FALSE]
  }
  
  mSigTools::write_exposure(
    exposure = final_sig_universe,
    file = file.path(synthetic_data_folder, dataset_name, "sig_universe_mask.csv")
  )
}

sort_non_msi_and_msi <-
  function(sample_names, sig_names) {
    
    # browser()
    
    out_matrix = matrix(0, nrow = length(sig_names), ncol = length(sample_names))
    colnames(out_matrix) = sample_names
    rownames(out_matrix) = sig_names
    
    # Check whether there are any MSI-H samples in spectra
    msi_samples <- grep(pattern = "MSI", x = sample_names)
    non_msi_sig_names <- setdiff(sig_names, get_msi_sig_names())
    
    if (length(msi_samples) > 0) {
      out_matrix[ , msi_samples] = 1
    }
    out_matrix[non_msi_sig_names, ] = 1 
    
    return(out_matrix)
  }

compute_sig_universe_mask("DBS")
compute_sig_universe_mask("SBS")
compute_sig_universe_mask("ID")