# Run this script with the top level directory as the working directory
source("analysis/code/analysis_utils.R")

output_home <- "analysis/raw_output/DBS/mp/syn"

time_used <- system.time({
  run_mp_syn(
    dataset_name = "DBS",
    output_home = output_home
  )
})

saveRDS(time_used, file = file.path(output_home, "time_used.Rds"))
