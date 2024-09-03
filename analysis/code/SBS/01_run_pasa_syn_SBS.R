# Run this script with the top level directory as the working directory
source("analysis/code/pasa_analysis.R")


run_pasa = function() {
  output_home <- "analysis/raw_output/SBS/pasa/syn"
  total_cores <- parallel::detectCores()
  mc_cores_per_sample <- min(5, total_cores)
  num_parallel_samples <- floor(total_cores / mc_cores_per_sample)
  
  time_used <- system.time({
    run_pasa_syn(
      dataset_name = "SBS",
      output_home = output_home,
      seed_in_use = 145879,
      num_parallel_samples = num_parallel_samples,
      mc_cores_per_sample = mc_cores_per_sample
    )
  })
  
  saveRDS(time_used, file = file.path(output_home, "time_used.Rds"))
}

run_pasa()

