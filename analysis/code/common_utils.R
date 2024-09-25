
get_tool_order <- function(assesment_by_sample, measure = mean) {
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


# These need to be the display tool names

global_tools_to_plot =
  c("PASA", "MuSiCal", "FitMS_01",  
    "SigPro", "MutPat", "YAPSA_03", "DeconSig_03",  "mutSig", "SigEstQP", "MSA_best", "sigLASSO", "sigLASSO w prior")


global_raw_tools_to_plot = 
  c("pasa", "musical", "fitms_0.010",  
    "sigpro", "mp", "yapsa_0.03", "deconstruct_0.03",  "mutsig", "sigest", 'msa_thresholdx100', "siglasso")

global_msa_index = which(global_raw_tools_to_plot == "msa_thresholdx100")
stopifnot(length(global_msa_index) == 1)

names(global_tools_to_plot) = global_raw_tools_to_plot

new_change_tool_names <- function(dt) {
  dt[Tool = global_tools_to_plot[Tool]]
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

plot_output_directory <- "output_for_paper/"
global_output_for_paper <- "output_for_paper/"

global_measures <- c("Combined", "one_minus_smd", "prec", "sens", "spec", "scaled_L2", "KL", "multiLLH")
global_ylabs = c(
  "Combined", "1 - scaled Manhattan distance",
  "Precision", "Recall (Sensitivity)", "Specificity", "1 - scaled L2", "log2(KL divergence + 1)", "Log likelihood")
names(global_ylabs) = global_measures
