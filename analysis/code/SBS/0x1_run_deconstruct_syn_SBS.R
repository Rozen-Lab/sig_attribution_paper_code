# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")

source("analysis/code/deconstruct_analysis.R")
source("analysis/code/generic_analysis.R")

run_deconstruct <- function() {
  
  more_args = list(signature_cutoff = 0.06) # The default for deconstructSigs
  
  output_home <- "analysis/raw_output/SBS/deconstruct/syn"
  time_used <- system.time({
    run_generic_syn(
      dataset_name = "SBS",
      output_home = output_home,
      attribute_function = call_deconstruct,
      more_args = more_args
    )
  })
  
  saveRDS(time_used, file = file.path(output_home, "time_used.Rds"))
}

run_deconstruct()
