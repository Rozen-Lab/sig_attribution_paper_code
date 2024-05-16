library(data.table)
library(writexl)
source("analysis/code/analysis_utils.R")

sbs_cpu_seconds_file <-
  "analysis/summary/SBS/syn/total_cpu_seconds.csv"
dbs_cpu_seconds_file <-
  "analysis/summary/DBS/syn/total_cpu_seconds.csv"
id_cpu_seconds_file <-
  "analysis/summary/ID/syn/total_cpu_seconds.csv"

sbs_df2 <- sbs_df <-
  data.table::fread(sbs_cpu_seconds_file)
sbs_df$total_cpu_seconds <- sbs_df$total_cpu_seconds / 60 / 60
sbs_df2$total_cpu_seconds <- sbs_df2$total_cpu_seconds / 900
colnames(sbs_df2)[2] <- colnames(sbs_df)[2] <- "SBS"

dbs_df2 <- dbs_df <-
  data.table::fread(dbs_cpu_seconds_file)
dbs_df$total_cpu_seconds <- dbs_df$total_cpu_seconds / 60 / 60
dbs_df2$total_cpu_seconds <- dbs_df2$total_cpu_seconds / 900
colnames(dbs_df2)[2] <- colnames(dbs_df)[2] <- "DBS"

id_df2 <- id_df <-
  data.table::fread(id_cpu_seconds_file)
id_df$total_cpu_seconds <- id_df$total_cpu_seconds / 60 / 60
id_df2$total_cpu_seconds <- id_df2$total_cpu_seconds / 900
colnames(id_df2)[2] <- colnames(id_df)[2] <- "ID"

df <- merge(x = sbs_df, y = dbs_df, by = "tool")
df2 <- merge(x = df, y = id_df, by = "tool")

cpu_seconds_per_sample <-
  merge(x = sbs_df2, y = dbs_df2, by = "tool")
cpu_seconds_per_sample2 <-
  merge(x = cpu_seconds_per_sample, y = id_df2, by = "tool")

sbs_file <-
  "analysis/summary/SBS/syn/assessment_each_sample.csv"
dt <- data.table::fread(sbs_file)
sbs_tool_order <- get_tool_order(dt)

df3 <- df2[match(sbs_tool_order, df2$tool), ]
colnames(df3)[1] <- "Tool"
df3 <- change_tool_names(df3)

cpu_seconds_per_sample3 <-
  cpu_seconds_per_sample2[match(sbs_tool_order, cpu_seconds_per_sample2$tool), ]
colnames(cpu_seconds_per_sample3)[1] <- "Tool"
cpu_seconds_per_sample3 <- change_tool_names(cpu_seconds_per_sample3)

#################################################################################
# Calculate elapsed time #
#################################################################################
sbs_elapsed_time_file <-
  "analysis/summary/SBS/syn/elapsed_time_seconds.csv"
dbs_elapsed_time_file <-
  "analysis/summary/DBS/syn/elapsed_time_seconds.csv"
id_elapsed_time_file <-
  "analysis/summary/ID/syn/elapsed_time_seconds.csv"

sbs_elapsed_df2 <- sbs_elapsed_df <-
  data.table::fread(sbs_elapsed_time_file)
sbs_elapsed_df2$elapsed_time_seconds <-
  sbs_elapsed_df2$elapsed_time_seconds / 900
colnames(sbs_elapsed_df2)[2] <- colnames(sbs_elapsed_df)[2] <- "SBS"

dbs_elapsed_df2 <- dbs_elapsed_df <-
  data.table::fread(dbs_elapsed_time_file)
dbs_elapsed_df2$elapsed_time_seconds <-
  dbs_elapsed_df2$elapsed_time_seconds / 900
colnames(dbs_elapsed_df2)[2] <- colnames(dbs_elapsed_df)[2] <- "DBS"

id_elapsed_df2 <- id_elapsed_df <-
  data.table::fread(id_elapsed_time_file)
id_elapsed_df2$elapsed_time_seconds <-
  id_elapsed_df2$elapsed_time_seconds / 900
colnames(id_elapsed_df2)[2] <- colnames(id_elapsed_df)[2] <- "ID"

elapsed_time_df <-
  merge(x = sbs_elapsed_df, y = dbs_elapsed_df, by = "tool")
elapsed_time_df2 <-
  merge(x = elapsed_time_df, y = id_elapsed_df, by = "tool")
elapsed_time_df3 <-
  elapsed_time_df2[match(sbs_tool_order, elapsed_time_df2$tool), ]
colnames(elapsed_time_df3)[1] <- "Tool"
elapsed_time_df3 <- change_tool_names(elapsed_time_df3)

elapsed_time_per_sample_df <-
  merge(x = sbs_elapsed_df2, y = dbs_elapsed_df2, by = "tool")
elapsed_time_per_sample_df2 <-
  merge(x = elapsed_time_per_sample_df, y = id_elapsed_df2, by = "tool")
elapsed_time_per_sample_df3 <-
  elapsed_time_per_sample_df2[match(
    sbs_tool_order,
    elapsed_time_per_sample_df2$tool
  ), ]
colnames(elapsed_time_per_sample_df3)[1] <- "Tool"
elapsed_time_per_sample_df3 <-
  change_tool_names(elapsed_time_per_sample_df3)

output_home <- "output_for_paper/"
writexl::write_xlsx(
  x = list(
    `Total CPU hours` = df3,
    `CPU seconds per sample` = cpu_seconds_per_sample3,
    `Elapsed seconds` = elapsed_time_df3,
    `Elapsed seconds per sample` = elapsed_time_per_sample_df3
  ),
  path = file.path(
    output_home,
    "sup_table_s8_syn_data_cpu_hours_elpased_time.xlsx"
  ),
  col_names = TRUE, format_headers = TRUE
)
