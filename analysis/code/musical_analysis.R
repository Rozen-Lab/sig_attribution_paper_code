source("analysis/code/generic_analysis.R")

library(reticulate)

# Run musical on one set of spectra with one set of signatures
call_musical <- function (spectra, 
                         signatures, 
                         more_args) {
  sigpro_path = file.path("analysis/raw_output", 
                          more_args$mutation_type, 
                          "musical/syn")
  msi_str     = ifelse(more_args$is_msi, "msi", "non_msi")
  
  # Need to put the .tsv files for one run of musical in a place from which python can read them
  input_dir  = file.path(sigpro_path, "input", more_args$cancer_type, msi_str)
  write_musical_input(spectra, signatures, input_dir)
  
  # expsoures = py_run_string('run_musical()')
  # or?
  browser()
  dictionary = py_run_string('exposures = run_musical()')
  
  return(dictionary)

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


NOT_NEEDED_EXCEPT_READ_TABLE_execute_sigpro_in_r <-
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

# This is a copy of write_sigpro_input
# Create tsv files for one run of musical
write_musical_input <- function(catalog, sig, output_dir) {
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


run_musical <- function(mut_type, more_args) {
  use_condaenv("musical2")
  source_python("analysis/code/run_musical.py")
  # browser()
  output_home <-
    file.path("analysis/raw_output", mut_type, "musical/syn")
  more_args$mutation_type = mut_type
  run_generic_syn(
    dataset_name = mut_type,
    output_home = output_home, # Duplicate info here and in call_sigpro
    attribute_function = call_musical,
    more_args = more_args
  )
}
