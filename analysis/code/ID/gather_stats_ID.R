source("analysis/code/analysis_utils.R")
library(parallel)

syn_exp_files <- list.files(
  path = "analysis/raw_output/ID/",
  full.names = TRUE, recursive = TRUE,
  pattern = "^inferred_exposures.csv"
)
tool_names <- basename(sub("/syn.*", "", syn_exp_files))

dataset <- "ID"
output_dir_syn <- "analysis/summary/ID/syn"

total_cores <- parallel::detectCores()
cores_to_use <- total_cores / 2

compare_syn_results(
  dataset = dataset,
  syn_exp_files = syn_exp_files,
  tool_names = tool_names,
  output_dir = output_dir_syn,
  mc_cores = cores_to_use
)

# foo2 <- read_csv("analysis/summary/DBS/syn/all_summary_stats.csv")