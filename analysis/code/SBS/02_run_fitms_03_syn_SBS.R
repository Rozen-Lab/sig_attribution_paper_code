# Run this script with the top level directory as the working directory
source("analysis/code/fitms_analysis.R")

run_fitms <- function() {
  
  output_home <- "analysis/raw_output/SBS/fitms_03/syn"
  time_used <- system.time({
    run_fitms_syn(
      dataset_name = "SBS",
      output_home = output_home,
      rare_sig_threshold = 0.03
    )
  })
  
  saveRDS(time_used, file = file.path(output_home, "time_used.Rds"))
}

run_fitms()
