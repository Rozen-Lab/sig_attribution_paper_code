source("analysis/code/analysis_utils.R")

output_home <- "analysis/raw_output/SBS/msa_opt/syn/output/"
mut_type <- "SBS"
summarize_msa_results(msa_output_dir = output_home, 
                      mut_type = mut_type)