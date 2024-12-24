global_random_seed = 145879


global_raw_tools_to_plot = 
  c("pasa", "musical", "fitms_0.010",  
    "sigpro", "mp", "yapsa", 
    "deconstruct",  "mutsig", "sigest",
    "siglasso", "sigfit", 
    "sigspack", "msa_default")

pretty_tool_names <- function(tool_names) {
  # browser()
  # These are the display tool names
  my_tools_to_plot =
    c("PASA", "MuSiCal", "FitMS",  
      "SigPro", "MutPat", "YAPSA", 
      "DeconSig",  "mutSig", "SigEstQP", 
      "sigLASSO",  "sigfit", 
      "SigsPack", "MSA")

  names(my_tools_to_plot) = global_raw_tools_to_plot

  retval = my_tools_to_plot[tool_names]
  # message("pretty in:  ", paste(tool_names, collapse = " "))
  # message("pretty out: ", paste(retval, collapse = " "))
  return(retval)
}


tools_to_plot_and_order = function(mutation_type) {
  filebase = paste0("summary_stats_", mutation_type, ".csv")
  dir = basename(getwd())
  if (dir == "output_for_paper") {
    file = filebase
  } else {
    file = file.path(global_output_for_paper, filebase)
  }
  ff = data.table::fread(file)
  ff= dplyr::filter(ff, Tool %in% global_raw_tools_to_plot)
  ff = dplyr::arrange(ff, desc(m.Combined))
  return(ff$Tool)
}


plot_output_directory <- "output_for_paper"
global_output_for_paper <- "output_for_paper"
global_summary_data_root = "output_for_paper"

global_measures <- c("Combined", "one_minus_smd", "prec", "sens", "spec", "scaled_L2", "KL")
global_ylabs = c(
  "Combined Score", "1 - scaled Manhattan distance",
  "Precision", "Recall (Sensitivity)", "Specificity", "1 - scaled L2 distance", "log2(KL divergence + 1)")
names(global_ylabs) = global_measures
