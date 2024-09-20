# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())

source("common_code/plot_functions.R")

data_home <- "analysis/summary/ID"
orig_indata <-
  data.table::fread(file.path(data_home, "assessment_each_sample_ID.csv"))
indata <- change_tool_names(orig_indata)

plot_objects <-
  boxplots_combined_cancer_types(
    assesment_by_sample = indata,
    tools_to_plot = global_tools_to_plot
  )

ggplot_to_pdf(
  plot_objects = plot_objects,
  file = "figure_5_ID_all_cancer_types.pdf",
  nrow = 4, ncol = 2,
  width = 8.2677, height = 11.6929, units = "in"
)
