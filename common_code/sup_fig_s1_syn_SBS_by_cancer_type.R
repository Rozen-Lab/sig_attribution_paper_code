# Run this script from the top level directory

source("analysis/code/analysis_utils.R")

data_home <- "analysis/summary/SBS/syn/"
indata <-
  data.table::fread(file.path(data_home, "assessment_each_sample.csv"))
indata <- change_tool_names(indata)

plot_objects <-
  boxplots_by_type(
    assesment_by_sample = indata
  )

output_home <- "output_for_paper/"

ggplot_to_pdf(
  plot_objects = plot_objects,
  file = file.path(
    output_home,
    "sup_fig_s1_syn_SBS_by_cancer_type.pdf"
  ),
  nrow = 5, ncol = 1,
  width = 8.2677, height = 11.6929, units = "in"
)
