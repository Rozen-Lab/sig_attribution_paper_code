source("analysis/code/analysis_utils.R")

total_cores <- parallel::detectCores()
cores_to_use <- total_cores / 2
output_home <- "analysis/raw_output/ID/"
total_cpu_seconds <-
  get_total_cpu_seconds(
    ouput_dir = output_home,
    mc_cores = cores_to_use
  )
elapsed_time <-
  get_elapsed_time(
    ouput_dir = output_home,
    mc_cores = cores_to_use
  )

output_dir_syn <- "analysis/summary/ID/syn"
summary_stats <-
  read.csv(file = file.path(output_dir_syn, "all_summary_stats.csv"))
summary_stats2 <- dplyr::arrange(summary_stats, desc(med.Combined))
tool_order <- summary_stats2$Tool
total_cpu_seconds2 <- total_cpu_seconds[tool_order]

df <- data.frame(
  tool = names(total_cpu_seconds2),
  total_cpu_seconds = total_cpu_seconds2
)
write.csv(
  x = df,
  file = file.path(output_dir_syn, "total_cpu_seconds.csv"),
  row.names = FALSE
)

elapsed_time2 <- elapsed_time[tool_order]
df2 <- data.frame(
  tool = names(elapsed_time2),
  elapsed_time_seconds = elapsed_time2
)
write.csv(
  x = df2,
  file = file.path(output_dir_syn, "elapsed_time_seconds.csv"),
  row.names = FALSE
)
