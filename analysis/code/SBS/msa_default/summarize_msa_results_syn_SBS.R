# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/msa_analysis.R")

output_home <- "analysis/raw_output/SBS/msa_default/syn/output/"
mut_type <- "SBS"
summarize_msa_results(msa_output_dir = output_home, 
                      mut_type = mut_type)
