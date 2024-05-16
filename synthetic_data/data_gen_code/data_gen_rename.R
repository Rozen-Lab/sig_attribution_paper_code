data_gen_rename <- 
  function(dataset, old_dataset_name, top_folder_name = "synthetic_data") {
  library(ICAMS)
  
  raw_syn_data_home <- file.path(top_folder_name, dataset, "intermed_results")

  final_syn_data_home <- file.path(top_folder_name, dataset)

  # standard base names of synthetic data files
  file_names <- c(
    "ground.truth.syn.catalog",
    "ground.truth.syn.exposures",
    "ground.truth.sig.universe",
    "ground.truth.sigs"
  )

  rename_from <- file.path(raw_syn_data_home, old_dataset_name)
  
  old_file_names <- list.files(path = rename_from, pattern = "\\.csv$")
  for (fn in file_names) {
    # Find the "old" file that matches the file name we want
    old_fn <- old_file_names[grepl(fn, old_file_names)]
    one_old_file <- file.path(rename_from, old_fn)
    one_new_file <- paste0(final_syn_data_home, .Platform$file.sep, fn, ".csv")
    message("rename ", one_old_file, " to ", one_new_file)
    rename_result <- 
      file.rename(from = one_old_file, to   = one_new_file)
    if (!all(rename_result)) {
      stop("Problem renaming ", one_old_file, " to ", one_new_file)
    }
  }
  
  unlink(rename_from, recursive = TRUE)
}
