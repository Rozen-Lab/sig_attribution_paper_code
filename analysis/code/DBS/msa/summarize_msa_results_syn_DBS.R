source("analysis/code/analysis_utils.R")

output_home <- "analysis/raw_output/DBS/msa/syn/output/"
mut_type <- "DBS"
summarize_msa_results(msa_output_dir = output_home, 
                      mut_type = mut_type)