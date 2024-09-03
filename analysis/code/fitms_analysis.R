source("analysis/code/analysis_utils.R")

# remotes::install_github("Nik-Zainal-Group/signature.tools.lib", ref = "v2.4.5")
# Note, as of 2024 09 03, the DESCRIPTION files still has Version: 2.4.4

library(signature.tools.lib)

run_fitms_msi <-
  function(spectra, msi_samples, common_sig, rare_sig, cancer_type) {
    all_retvals <- list()
    common_sig_names <- colnames(common_sig)
    rare_sig_names <- colnames(rare_sig)
    
    if (length(msi_samples) == 0) {
      retval <-
        signature.tools.lib::FitMS(
          catalogues = spectra,
          commonSignatures = common_sig,
          rareSignatures = rare_sig
        )
      all_retvals[[cancer_type]] <- retval
    } else {
      spectra_msi <- spectra[, msi_samples, drop = FALSE]
      common_sig_msi <- common_sig
      rare_sig_msi <- rare_sig
      
      non_msi_samples <- setdiff(colnames(spectra), msi_samples)
      spectra_non_msi <- spectra[, non_msi_samples, drop = FALSE]
      msi_sig_names <- get_msi_sig_names()
      common_sig_names_non_msi <- setdiff(common_sig_names, msi_sig_names)
      rare_sig_names_non_msi <- setdiff(rare_sig_names, msi_sig_names)
      common_sig_non_msi <-
        common_sig[, common_sig_names_non_msi, drop = FALSE]
      rare_sig_non_msi <- rare_sig[, rare_sig_names_non_msi, drop = FALSE]
      
      retval_msi <-
        signature.tools.lib::FitMS(
          catalogues = spectra_msi,
          commonSignatures = common_sig_msi,
          rareSignatures = rare_sig_msi
        )
      all_retvals[[paste0(cancer_type, "_msi")]] <- retval_msi
      
      retval_non_msi <-
        signature.tools.lib::FitMS(
          catalogues = spectra_non_msi,
          commonSignatures = common_sig_non_msi,
          rareSignatures = rare_sig_non_msi
        )
      all_retvals[[paste0(cancer_type, "_non_msi")]] <- retval_non_msi
    }
    return(all_retvals)
  }

run_fitms_syn <- function(dataset_name, output_home,
                          data_top_folder_name = "synthetic_data",
                          rare_sig_threshold = 0.01,
                          cancer_types = NULL) {
  message("Start running the job")
  message("rare_sig_threshold = ", rare_sig_threshold)
  
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
    sig_props <- sig_universe[[cancer_type]]
    common_sig_names <- names(sig_props[sig_props >= rare_sig_threshold])
    rare_sig_names <- names(sig_props[sig_props < rare_sig_threshold])
    
    common_sig <- gt_sig[, common_sig_names, drop = FALSE]
    rare_sig <- gt_sig[, rare_sig_names, drop = FALSE]
    
    # Check whether there are any MSI-H samples in spectra
    msi_samples <- grep(pattern = "MSI", x = colnames(spectra), value = TRUE)
    
    time_used <- system.time({
      retval <- run_fitms_msi(
        spectra = spectra, msi_samples = msi_samples,
        common_sig = common_sig, rare_sig = rare_sig,
        cancer_type = cancer_type
      )
    }) # system.time
    
    time_by_cancer_type[[cancer_type]] <- time_used    
    all_retvals <- c(all_retvals, retval)
    
  } # for (cancer_type in cancer_types...
  
  all_exposures <- lapply(all_retvals, FUN = function(retval) {
    exposure <- t(retval$exposures)
    # Remove the row containing unassigned mutations
    exposure2 <- exposure[rownames(exposure) != "unassigned", , drop = FALSE]
    return(exposure2)
  })
  
  saveRDS(time_by_cancer_type,
          file.path(output_home, "time_by_cancer_type.Rds"))
  
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

