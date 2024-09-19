

run_generic_msi <-
  function(spectra, sig, cancer_type, more_args, attribute_function) {
    all_retvals <- list()
    msi_sig_names <- get_msi_sig_names()
    
    # Check whether there are any MSI-H samples in spectra
    msi_samples <- grep(pattern = "MSI", x = colnames(spectra), value = TRUE)

    if (length(msi_samples) > 0) {      
      spectra_msi <- spectra[, msi_samples, drop = FALSE]

      non_msi_samples <- setdiff(colnames(spectra), msi_samples)
      spectra_non_msi <- spectra[, non_msi_samples, drop = FALSE]
      non_msi_sig_names <- setdiff(colnames(sig), msi_sig_names)
      
      more_args$is_msi <- TRUE

      retval_msi <-
        attribute_function(
          spectra = spectra_msi,
          signatures = sig, # All signatures, including MSI signatures
          more_args
        )
      
      all_retvals[[paste0(cancer_type, "_msi")]] <- retval_msi
      retval_key <- paste0(cancer_type, "_non_msi")
      
    } else {
      non_msi_sig_names <- colnames(sig)
      non_msi_samples <- colnames(spectra)
      retval_key <- cancer_type
    }
    
    more_args$is_msi <- FALSE
    
    retval_non_msi <-
      attribute_function(
        spectra = spectra[ , non_msi_samples, drop = FALSE],
        signatures = sig[, non_msi_sig_names, drop = FALSE],
        more_args
      )
    
    all_retvals[[retval_key]] <- retval_non_msi

    return(all_retvals)
  }



run_generic_syn <- function(dataset_name, output_home,
                                data_top_folder_name = "synthetic_data",
                                cancer_types = NULL,
                                attribute_function,
                                more_args) {
  message("Start running the job")

  if (!dir.exists(output_home)) {
    dir.create(path = output_home, recursive = TRUE)
  }
  
  all_inputs <- get_all_input( # defined in analysis_utils.R
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
  
  outer_time_used = system.time({
    for (cancer_type in cancer_types) { 
      message(cancer_type)
      spectra <- all_spectra[[cancer_type]]
      sig_names <- names(sig_universe[[cancer_type]])
      ####
      more_args$sig_props <- sig_universe[[cancer_type]]
      more_args$cancer_type <- cancer_type
      ####
      sig <- gt_sig[, sig_names, drop = FALSE]
      output_dir <- file.path(output_home, cancer_type)
      
      inner_time_used <- system.time({
        
        retval <-
          run_generic_msi(
            spectra = spectra,
            sig = sig,
            cancer_type = cancer_type,
            attribute_function = attribute_function,
            more_args
          )
        # retval is a list of matrices of with rows being signatures and columns being samples
      }) # system.time
      
      time_by_cancer_type[[cancer_type]] <- inner_time_used
      all_retvals <- c(all_retvals, retval)
      
    } # for (cancer_type in cancer_types...
  }) # outer_time_used
  
  file = file.path(output_home, "time_used.Rds")
  message("Saving ", file)
  saveRDS(outer_time_used, file = file.path(output_home, "time_used.Rds"))
  
  file = file.path(output_home, "time_by_cancer_type.Rds")
  message("Saving ", file)
  saveRDS(time_by_cancer_type, file)
  
  final_exposure <- mSigAct:::MergeListOfExposures(all_retvals)
  
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
  
  file = file.path(output_home, "all_retvals.Rds")
  message("Saving ", file)
  saveRDS(all_retvals, file = file)
  message("Finished running all the jobs")
}


get_all_input <- function(dataset_name,
                          data_top_folder_name = "synthetic_data") {
  data_home <- file.path(data_top_folder_name, dataset_name)
  
  spectra <- ICAMS::ReadCatalog(
    file = file.path(data_home, "ground.truth.syn.catalog.csv")
  )
  
  ground_truth_sigs <- ICAMS::ReadCatalog(
    file.path(data_home, "ground.truth.sigs.csv"),
    catalog.type = "counts.signature"
  )
  
  sig_universe_info <- data.table::fread(
    file.path(data_home, "ground.truth.sig.universe.csv"),
    fill = TRUE
  )
  
  spectra_list <- PCAWG7::SplitPCAWGMatrixByTumorType(spectra)
  
  cancer_types <- names(spectra_list)
  
  signature_universes <-
    sapply(cancer_types, FUN = function(cancer_type) {
      one_type_info <- sig_universe_info[V1 == cancer_type, -1]
      sig_props <- unlist(one_type_info[2, ])
      sig_props <- as.numeric(sig_props[sig_props != ""])
      
      sig_names <- unlist(one_type_info[1, ])
      sig_names <- sig_names[sig_names != ""]
      names(sig_props) <- sig_names
      return(sig_props)
    })
  
  return(list(
    spectra_list = spectra_list,
    ground_truth_sigs = ground_truth_sigs,
    signature_universes = signature_universes,
    cancer_types = cancer_types
  ))
}

get_msi_sig_names <- function() {
  sig_names <-
    c("SBS6", "SBS14", "SBS15", "SBS20", "SBS21", "SBS26", "SBS44")
  return(sig_names)
}


