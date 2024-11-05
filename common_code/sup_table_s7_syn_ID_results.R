library(writexl)
library(data.table)
source("analysis/code/analysis_utils.R")

file <-
  "analysis/summary/ID/syn/assessment_each_sample.csv"
dt <- data.table::fread(file)
stats_table <- create_median_iqr_table(dt)

output_home <- "output_for_paper/"
writexl::write_xlsx(
  x = list(`Table S7` = stats_table),
  path = file.path(output_home, "sup_table_s7_syn_ID_results.xlsx"),
  col_names = TRUE, format_headers = TRUE
)
