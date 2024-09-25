# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())

source("analysis/code/plot_functions.R")

plot_by_cancer_type_sup_fig("ID", sup_fig_num = "3")