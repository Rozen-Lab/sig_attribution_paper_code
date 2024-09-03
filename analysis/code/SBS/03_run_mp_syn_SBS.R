# Run this script with the top level directory as the working directory
source("analysis/code/mp_analysis.R")

run_mp = function() {
  output_home <- "analysis/raw_output/SBS/mp/syn"
  
  time_used <- system.time({
    run_mp_syn(
      dataset_name = "SBS",
      output_home = output_home
    )
  })
  
  saveRDS(time_used, file = file.path(output_home, "time_used.Rds"))
}

run_mp()
