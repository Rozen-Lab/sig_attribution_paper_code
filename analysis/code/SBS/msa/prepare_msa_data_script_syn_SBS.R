source("analysis/code/analysis_utils.R")

dataset_name <- "SBS"
input_dir <- "analysis/raw_output/SBS/msa/syn/input"
data_identifier <- "syn_data"
output_home <- "analysis/raw_output/SBS/msa/syn/output"
bash_script_home <- "analysis/code/SBS/msa/"
msa_nextflow_file <- "common_code/MSA/run_auto_optimised_analysis.nf"
project_dir <- getwd()

prepare_msa_data_script_syn(
  dataset_name = dataset_name,
  input_dir = input_dir,
  data_identifier = data_identifier,
  output_home = output_home,
  bash_script_home = bash_script_home,
  msa_nextflow_file = msa_nextflow_file,
  project_dir = project_dir
)
