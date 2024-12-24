source("analysis/code/generic_analysis.R")

# This code requires conda and the musical2 conda environment; see
# file conda_env_for_musical.yml in this directory
# You will need to edit the line in the .yml file that points
# to the environment (the last line in the file)

# Run musical on one set of spectra with one set of signatures
call_musical <- function (spectra, 
                         signatures, 
                         more_args) {
  

  musical_path = file.path("analysis/raw_output", 
                          more_args$mutation_type, 
                          "musical/syn")
  msi_str     = ifelse(more_args$is_msi, "msi", "non_msi")
  
  # Need to put the .tsv files for one run of musical in a place from which python can read them
  input_dir  = file.path(musical_path, "input", more_args$cancer_type, msi_str)
  output_dir  = file.path(musical_path, "output", more_args$cancer_type, msi_str)
  write_musical_input(spectra, signatures, input_dir)
  print(input_dir)
  print(output_dir)
  stopifnot(file.exists("analysis/code/run_musical.py"))
  # browser()
  pyreturn = system2("conda",  c("run", "-n",  "musical2", "python", "analysis/code/run_musical.py", input_dir, output_dir, "20240916"))
  if (pyreturn != 0) {
    stop("There was an error running the python code")
  }
  exposure = read.csv(file.path(output_dir, "exposures.csv"), row.names = 1, check.names = FALSE)
  return(exposure)
                            
} 

# For musical the code needs to create and save inputs to
# analysis/raw_output/<mut_type>/musical/syn/input/<cancertype>/{msi,non-msi}/{catalog.tsv, sigs.tsv}
#
# outputs go to
# analysis/raw_output/<mut_type>/musical/syn/output/, 
# with musical outputs in folders named by cancer type.
# each cancer-type specific folder contains non-msi and
# possible msi folder, which contains musical output


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
  # browser()
  testconda = system2("conda",  c("run", "-n",  "musical2", "echo", "ok"))
  if (testconda != 0) {
    stop("There was an error activating the musical2 conda enviroment\n",
         "Please create the musical2 conda environment from conda_env_for_musical.yml")
  }
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
