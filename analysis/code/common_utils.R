global_random_seed = 145879

get_tool_order <- function(assesment_by_sample, measure = mean) {
  # browser()
  measure_combined <-
    aggregate(assesment_by_sample$Combined,
              list(assesment_by_sample$Tool),
              FUN = measure
    )
  tool_order <-
    measure_combined[order(measure_combined$x, decreasing = TRUE), ]$Group.1
  return(tool_order)
}

# Not sure if this is used
if (FALSE) {
all_tools = 
  c("PASA", "musical", "fitms_0.010", "fitms_0.030", "fitms_0.060",  "fitms_0.005", 
    "SigPro", "MP", "yapsa_0.03", "deconstruct_0.03",  "yapsa_0.01", "deconstruct_0.06",
    "yapsa_0.06", "mutsig", "yapsa_0.00",  "deconstruct_0.10", "deconstruct_0.00", "yapsa_0.10",
    "QP_null",  "sigest", "msa", "msa_unpruned_default")
}



pretty_tool_names <- function(tool_names) {
  
  # These are the display tool names
  my_tools_to_plot =
    c("PASA", "MuSiCal", "FitMS_01",  
      "SigPro", "MutPat", "YAPSA_03", "YAPSA_06", 
      "DeconSig_03",  "mutSig", "SigEstQP", 
      "MSA_best", "sigLASSO",  "sigfit", "SigsPack" )
  my_raw_tools_to_plot = 
    c("pasa", "musical", "fitms_0.010",  
      "sigpro", "mp", "yapsa_0.03", "yapsa_0.06", 
      "deconstruct_0.03",  "mutsig", "sigest",
      'msa_thresholdx100', "siglasso", "sigfit", "sigspack")
  names(my_tools_to_plot) = my_raw_tools_to_plot

  retval = my_tools_to_plot[tool_names]
  return(retval)
}

change_tool_names <- function(dt) {
  dt[Tool == "pasa", Tool := "PASA"]
  dt[Tool == "fitms_0.010", Tool := "FitMS_01"]
  dt[Tool == "sigpro", Tool := "SigPro"]
  dt[Tool == "mp", Tool := "MutPat"]
  dt[Tool == "msa", Tool := "MSA"]
  dt[Tool == "msa_opt", Tool := "MSA_opt"]
  dt[Tool == "yapsa_0.03", Tool := "YAPSA_03"]
  dt[Tool == "deconstruct_0.03", Tool := "DeconSig_03"]
  dt[Tool == "mutsig", Tool := "mutSig"]
  dt[Tool == "sigest", Tool := "SigEstQP"]
  dt[Tool == "musical", Tool := "MuSiCal"]
  dt[Tool == "siglasso", Tool := "sigLASSO"]
  dt[Tool == "siglasso_wprior", Tool := "sigLASSO w prior"]
  dt[Tool == "msa_default_unpruned", Tool := "MSA"]
  return(dt)
}

raw_tools_to_plot = function(mutation_type) {
  if (mutation_type == "SBS") {
    return(
      c("pasa",
        "musical",
        "fitms_0.010",  
        "sigpro",
        "mp",
        "yapsa_0.03", 
        "deconstruct_0.03",
        "mutsig",
        "sigest",
        'msa_thresholdx100',
        "siglasso", 
        "sigfit",
        "sigspack"))
  } else if (mutation_type == "DBS") {
    return(
      c("pasa",
        "fitms_0.010",
        "musical",
        "mp",
        "msa_thresholdx100",
        "yapsa_0.06", # Different from SBS
        "deconstruct_0.03",
        "mutsig",
        "sigpro",
        "sigfit",
        "sigspack",
        "sigest"
      )
    )
  } else if (mutation_type == "ID") {
    return(
      c("pasa",
        "fitms_0.010",
        "musical",
        "mp",
        "msa_thresholdx100",
        "yapsa_0.06", # Same as DBS, different from SBS
        "deconstruct_0.03",
        "mutsig",
        "sigpro",
        "sigfit",
        "sigspack",
        "sigest"
      )
    )
  } else {
    stop("Unkown mutation_type: ", mutation_type)
  }
  
}


plot_output_directory <- "output_for_paper/"
global_output_for_paper <- "output_for_paper/"

global_measures <- c("Combined", "one_minus_smd", "prec", "sens", "spec", "scaled_L2", "KL", "multiLLH")
global_ylabs = c(
  "Combined", "1 - scaled Manhattan distance",
  "Precision", "Recall (Sensitivity)", "Specificity", "1 - scaled L2", "log2(KL divergence + 1)", "Log likelihood")
names(global_ylabs) = global_measures
