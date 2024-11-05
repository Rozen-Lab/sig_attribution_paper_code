source("analysis/code/generic_analysis.R")

call_pasa <- function (spectra, 
                       signatures, 
                       more_args
) {
  # browser()
  retval <- mSigAct::PresenceAttributeSigActivity(
    spectra                   = spectra,
    sigs                      = signatures,
    output.dir                = more_args$output_dir,
    num.parallel.samples      = more_args$num_parallel_sample,
    mc.cores.per.sample       = more_args$mc_cores_per_sample,
    seed                      = more_args$seed_in_use,
    save.files                = FALSE
  )
  pasa_report_error_msg(retval = retval, cancer_type = cancer_type)
  # browser()
  exposure_matrix = retval$proposed.assignment
  return(exposure_matrix)
}


run_pasa <- function(mut_type, more_args) {
  output_home <-
    file.path("analysis/raw_output", mut_type, "pasa/syn")
  message("PASA, writing to ", output_home)
  more_args$output_home = output_home
  
  run_generic_syn(
    dataset_name = mut_type,
    output_home = output_home,
    attribute_function = call_pasa,
    more_args = more_args
  )
}


pasa_report_error_msg <- function(retval, cancer_type) {
  errs <- retval$error.messages
  err_index <- which(errs != "")
  if (length(err_index) > 0) {
    message("Got errors for cancer type ", cancer_type)
    message(paste(names(errs)[err_index], collapse = " "))
    message(paste(errs[err_index], collapse = " "))
  }
}
