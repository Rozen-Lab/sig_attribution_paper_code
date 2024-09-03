source("analysis/code/analysis_utils.R")

pasa_report_error_msg <- function(retval, cancer_type) {
  errs <- retval$error.messages
  err_index <- which(errs != "")
  if (length(err_index) > 0) {
    message("Got errors for cancer type ", cancer_type)
    message(paste(names(errs)[err_index], collapse = " "))
    message(paste(errs[err_index], collapse = " "))
  }
}

execute_pasa <-
  function(spectra, sig, output_dir, num_parallel_samples,
           mc_cores_per_sample, seed_in_use, cancer_type) {
    retval <- mSigAct::PresenceAttributeSigActivity(
      spectra                   = spectra,
      sigs                      = sig,
      output.dir                = output_dir,
      num.parallel.samples      = num_parallel_samples,
      mc.cores.per.sample       = mc_cores_per_sample,
      seed                      = seed_in_use,
      save.files                = FALSE
    )
    pasa_report_error_msg(retval = retval, cancer_type = cancer_type)
    return(retval)
  }

run_pasa_msi <-
  function(spectra,
           sig, 
           output_dir,
           num_parallel_samples, 
           mc_cores_per_sample, 
           seed_in_use,
           cancer_type) {
    
    # Check whether there are any MSI-H samples in spectra
    msi_samples <- grep(pattern = "MSI", x = colnames(spectra), value = TRUE)
    msi_sig <- get_msi_sig_names()
    all_retvals <- list()
    if (length(msi_samples) == 0) {
      retval <-
        execute_pasa(
          spectra = spectra, sig = sig, output_dir = output_dir,
          num_parallel_samples = num_parallel_samples,
          mc_cores_per_sample = mc_cores_per_sample,
          seed_in_use = seed_in_use,
          cancer_type = cancer_type
        )
      all_retvals[[cancer_type]] <- retval
    } else {
      spectra_msi <- spectra[, msi_samples, drop = FALSE]
      sig_msi <- sig
      output_dir_msi <- file.path(output_dir, "msi")

      non_msi_samples <- setdiff(colnames(spectra), msi_samples)
      spectra_non_msi <- spectra[, non_msi_samples, drop = FALSE]
      non_msi_sig_names <- setdiff(colnames(sig), msi_sig)
      sig_non_msi <- sig[, non_msi_sig_names, drop = FALSE]
      output_dir_non_msi <- file.path(output_dir, "non_msi")

      retval_msi <-
        execute_pasa(
          spectra = spectra_msi, sig = sig_msi,
          output_dir = output_dir_msi,
          num_parallel_samples = num_parallel_samples,
          mc_cores_per_sample = mc_cores_per_sample,
          seed_in_use = seed_in_use,
          cancer_type = cancer_type
        )
      all_retvals[[paste0(cancer_type, "_msi")]] <- retval_msi

      retval_non_msi <-
        execute_pasa(
          spectra = spectra_non_msi, sig = sig_non_msi,
          output_dir = output_dir_non_msi,
          num_parallel_samples = num_parallel_samples,
          mc_cores_per_sample = mc_cores_per_sample,
          seed_in_use = seed_in_use,
          cancer_type = cancer_type
        )
      all_retvals[[paste0(cancer_type, "_non_msi")]] <- retval_non_msi
    }
    return(all_retvals)
  }

run_pasa_syn <- function(dataset_name, output_home,
                         seed_in_use,
                         data_top_folder_name = "synthetic_data",
                         num_parallel_samples = 1,
                         mc_cores_per_sample = 5,
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

    time_used <- system.time({
      
      retval <-
        run_pasa_msi(
          spectra = spectra,
          sig = sig, 
          output_dir = output_dir,
          num_parallel_samples = num_parallel_samples,
          mc_cores_per_sample = mc_cores_per_sample,
          seed_in_use = seed_in_use,
          cancer_type = cancer_type
        )
      
    }) # system.time
    
    time_by_cancer_type[[cancer_type]] <- time_used
    all_retvals <- c(all_retvals, retval)
    
  } # for (cancer_type in cancer_types)
  
  saveRDS(time_by_cancer_type,
          file.path(output_home, "time_by_cancer_type.Rds"))
  
  all_exposures <- lapply(all_retvals, FUN = function(retval) {
    return(retval$proposed.assignment)
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

