source("analysis/code/analysis_utils.R")

library(MutationalPatterns)

run_mp_msi <-
  function(spectra, msi_samples, sig, cancer_type) {
    all_retvals <- list()
    msi_sig_names <- get_msi_sig_names()
    
    if (length(msi_samples) == 0) {
      retval <-
        MutationalPatterns::fit_to_signatures_strict(
          mut_matrix = spectra,
          signatures = sig
        )
      all_retvals[[cancer_type]] <- retval
    } else {
      spectra_msi <- spectra[, msi_samples, drop = FALSE]
      sig_msi <- sig
      
      non_msi_samples <- setdiff(colnames(spectra), msi_samples)
      spectra_non_msi <- spectra[, non_msi_samples, drop = FALSE]
      non_msi_sig_names <- setdiff(colnames(sig), msi_sig_names)
      sig_non_msi <- sig[, non_msi_sig_names, drop = FALSE]
      
      retval_msi <-
        MutationalPatterns::fit_to_signatures_strict(
          mut_matrix = spectra_msi,
          signatures = sig_msi
        )
      all_retvals[[paste0(cancer_type, "_msi")]] <- retval_msi
      
      retval_non_msi <-
        MutationalPatterns::fit_to_signatures_strict(
          mut_matrix = spectra_non_msi,
          signatures = sig_non_msi
        )
      all_retvals[[paste0(cancer_type, "_non_msi")]] <- retval_non_msi
    }
    return(all_retvals)
  }

run_mp_syn <- function(dataset_name, output_home,
                       data_top_folder_name = "synthetic_data",
                       cancer_types = NULL) {
  message("Start running the job")
  
  if (!dir.exists(output_home)) {
    dir.create(path = output_home, recursive = TRUE)
  }
  
  all_inputs <- get_all_input(
    dataset_name = dataset_name,
    data_top_folder_name = data_top_folder_name
  )
  all_spectra <- all_inputs$spectra_list
  gt_sig <- all_inputs$ground_truth_sigs
  sig_universe <- all_inputs$signature_universes
  
  if (is.null(cancer_types)) {
    cancer_types <- all_inputs$cancer_types
  }
  
  all_retvals <- list()
  time_by_cancer_type <- list()
  for (cancer_type in cancer_types) {
    spectra <- all_spectra[[cancer_type]]
    sig_names <- names(sig_universe[[cancer_type]])
    sig <- gt_sig[, sig_names, drop = FALSE]
    output_dir <- file.path(output_home, cancer_type)
    
    # Check whether there are any MSI-H samples in spectra
    msi_samples <- grep(pattern = "MSI", x = colnames(spectra), value = TRUE)
    
    time_used <- system.time({
      
      retval <-
        run_mp_msi(
          spectra = spectra,
          msi_samples = msi_samples,
          sig = sig,
          cancer_type = cancer_type
        )
    }) # system.time
    
    time_by_cancer_type[[cancer_type]] <- time_used
    all_retvals <- c(all_retvals, retval)
    
  } # for (cancer_type in cancer_types...
  
  saveRDS(time_by_cancer_type,
          file.path(output_home, "time_by_cancer_type.Rds"))
  
  all_exposures <- lapply(all_retvals, FUN = function(retval) {
    return(mSigAct:::RemoveZeroActivitySig(retval$fit_res$contribution))
  })
  
  final_exposure <- mSigAct:::MergeListOfExposures(all_exposures)
  
  original_spectra_file <-
    file.path(data_top_folder_name, dataset_name, "ground.truth.syn.catalog.csv")
  
  original_spectra <- ICAMS::ReadCatalog(original_spectra_file)
  
  if (setequal(colnames(original_spectra), colnames(final_exposure))) {
    final_exposure <-
      final_exposure[, colnames(original_spectra), drop = FALSE]
  }
  
  mSigTools::write_exposure(
    exposure = final_exposure,
    file = file.path(output_home, "inferred_exposures.csv")
  )
  
  saveRDS(all_retvals, file = file.path(output_home, "all_retvals.Rds"))
  message("Finished running all the jobs")
}

