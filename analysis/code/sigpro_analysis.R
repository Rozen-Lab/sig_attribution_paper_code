source("analysis/code/generic_analysis.R")
source("analysis/code/common_utils.R")

call_sigpro <- function (spectra, 
                         signatures, 
                         more_args) {
  sigpro_path = file.path("analysis/raw_output", 
                          more_args$mutation_type, 
                          "sigpro/syn")
  msi_str     = ifelse(more_args$is_msi, "msi", "non_msi")
  
  # Need to put the files in sigpro input dir!!!  
  sigpro_input_dir       = file.path(sigpro_path, "input", more_args$cancer_type, msi_str)
  sigpro_output_dir      = file.path(sigpro_path, "output", more_args$cancer_type, msi_str)
  
  write_sigpro_input(spectra, signatures, sigpro_input_dir)
  
  sigpro_exposure <-
    execute_sigpro_in_r(
      python_bin      = more_args$python_bin,
      run_sigpro_file = "analysis/code/run_sigpro.py",
      input_dir       = sigpro_input_dir,
      output_dir      = sigpro_output_dir,
      seed_in_use     = more_args$seed_in_use,
      context_type    = more_args$context_type
    )
  return(sigpro_exposure)
} 


# For sigpro we need to create and save inputs to
# analysis/raw_output/<mut_type>/sigpro/syn/input/<cancertype>/{msi,non-msi}/{catalog.tsv, sigs.tsv}
#
# outputs go to
# analysis/raw_output/<mut_type>/sigpro/syn/output/, 
# with sipro specific outputs in folders named by cancer type and
# regular outputs in at the top level of that folder
# each cancer-type specific (sigpro-speical) folder contains non-msi and
# possible msi folder, whicih contains sigpro output


execute_sigpro_in_r <-
  function(python_bin, run_sigpro_file, input_dir,
           output_dir, seed_in_use, context_type) {
    py_args <-
      c(run_sigpro_file, input_dir, output_dir, seed_in_use, context_type)
    system_ret <- system2(python_bin, args = py_args)
    if (system_ret != 0) {
      msg <- paste(py_args, collapse = " ")
      stop("Could not run this command:\n", python_bin, " ", msg)
    }
    
    sigpro_output <-
      file.path(
        output_dir,
        "Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
      )
    sigpro_exposure <-
      t(read.table(sigpro_output, header = TRUE, sep = "\t", row.names = 1))
    return(sigpro_exposure)
  }


write_sigpro_input <- function(catalog, sig, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  ICAMS:::ConvertCatalogToSigProfilerFormat(
    input.catalog = catalog,
    file = file.path(output_dir, "catalog.tsv")
  )
  ICAMS:::ConvertCatalogToSigProfilerFormat(
    input.catalog = sig,
    file = file.path(output_dir, "sigs.tsv")
  )
}


run_sigpro <- function(mut_type, more_args) {
  # browser()
  output_home <-
    file.path("analysis/raw_output", mut_type, "sigpro/syn")
  more_args$mutation_type = mut_type
  run_generic_syn(
    dataset_name = mut_type,
    output_home = output_home, # Duplicate info here and in call_sigpro
    attribute_function = call_sigpro,
    more_args = more_args
  )
}

