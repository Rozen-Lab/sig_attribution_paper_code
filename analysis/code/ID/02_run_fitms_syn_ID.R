# Run this script with the top level directory as the working directory
source("analysis/code/analysis_utils.R")

output_home <- "analysis/raw_output/ID/fitms/syn"

time_used <- system.time({
  run_fitms_syn(
    dataset_name = "ID",
    output_home = output_home
  )
})

saveRDS(time_used, file = file.path(output_home, "time_used.Rds"))
