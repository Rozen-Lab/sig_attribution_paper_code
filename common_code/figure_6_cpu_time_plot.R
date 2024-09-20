# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())

source("common_code/plot_functions.R")

sbs_cpu_seconds_file <-
  "analysis/summary/SBS/total_cpu_seconds_SBS.csv"
sbs_cpu_plot_object <-
  cpu_barplot(sbs_cpu_seconds_file,
              main = "Synthetic SBS",
              ylim = c(0, 450)
  )

dbs_cpu_seconds_file <-
  "analysis/summary/DBS/total_cpu_seconds_DBS.csv"
dbs_cpu_plot_object <-
  cpu_barplot(dbs_cpu_seconds_file,
              main = "Synthetic DBS",
              ylim = c(0, 450)
  )

id_cpu_seconds_file <-
  "analysis/summary/ID/total_cpu_seconds_ID.csv"
# browser()
id_cpu_plot_object <-
  cpu_barplot(id_cpu_seconds_file,
              main = "Synthetic ID",
              ylim = c(0, 450)
  )

output_home <- "output_for_paper/"

ggplot_to_pdf(
  plot_objects = c(
    list(sbs_cpu_plot_object),
    list(dbs_cpu_plot_object),
    list(id_cpu_plot_object)
  ),
  file = "figure_6_cpu_time_three_types.pdf",
  nrow = 3, ncol = 3,
  width = 8.2677, height = 11.6929, units = "in"
)
