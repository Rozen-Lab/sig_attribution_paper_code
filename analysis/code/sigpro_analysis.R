source("analysis/code/generic_analysis.R")

call_sigpro <- function (spectra, 
                         signatures, 
                         more_args
) {
  # browser()
    retval = run_sigpro_syn(
      dataset_name = "DBS",
      python_bin = "/home/e0012078/software/miniconda3/bin/python",
      run_sigpro_file = "analysis/code/run_sigpro.py",
      seed_in_use     = more_args$seed_in_use, #145879,
      context_type    = more_args$context_type, # "DINUC",
      input_root = "analysis/raw_output/DBS/sigpro/syn/input",
      output_root = "analysis/raw_output/DBS/sigpro/syn/output"
    )
  })
  return(retval)
}

run_sigpro <- function(mut_type, more_args) {
  out_fragment = paste0("fitms_", formatC(rare_sig_threshold, digits = 2, format = "f"))
  more_args$ <- list(rare_sig_threshold = rare_sig_thres
  output_home <-
    file.path("analysis/raw_output", mut_type, "sigpro/syn")
  time_used <- system.time({
    run_generic_syn(
      dataset_name = mut_type,
      output_home = output_home,
      attribute_function = call_sigpro,
      more_args = more_args
    )
  })
  
  saveRDS(time_used, file = file.path(output_home, "time_used.Rds"))
}

run_sigpro <- function(mut_type, more_args) {
  output_home <-
    file.path("analysis/raw_output", mut_type, "sigpro/syn")
  message("sigpro, writing to ", output_home)
  more_args$output_home = output_home
  time_used <- system.time({
    run_generic_syn(
      dataset_name = mut_type,
      output_home = output_home,
      attribute_function = call_sigpro,
      more_args = more_args
    )
  })
  
  saveRDS(time_used, file = file.path(output_home, "time_used.Rds"))
}

