# Run this script with the top level directory as the working directory
# For each mutation type (SBS, DBS, ID) this script creates one main text 
# figure and one supplementary figure
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/plot_functions.R")
plot_all_cancer_types_merged("SBS", "3", sup_fig_num = "1")
plot_all_cancer_types_merged("DBS", "4", sup_fig_num = "4")
plot_all_cancer_types_merged("ID",  "5", sup_fig_num = "6")
