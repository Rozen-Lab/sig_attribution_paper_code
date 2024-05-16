# Run this script from the top level directory

source("analysis/code/analysis_utils.R")

data_home <- "analysis/summary/ID/syn/"
indata <-
  data.table::fread(file.path(data_home, "assessment_each_sample.csv"))
indata <- change_tool_names(indata)

plot_objects <-
  boxplots_all_types(
    assesment_by_sample = indata
  )

output_home <- "output_for_paper/"

ggplot_to_pdf(
  plot_objects = plot_objects,
  file = file.path(
    output_home,
    "figure_5_syn_ID_all_cancer_types.pdf"
  ),
  nrow = 4, ncol = 2,
  width = 8.2677, height = 11.6929, units = "in"
)
