# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/generic_analysis.R")
source("analysis/code/msa_analysis.R")

dataset_name <- "ID"
input_dir <- "analysis/raw_output/ID/msa_default/syn/input"
data_identifier <- "syn_data"
output_home <- "analysis/raw_output/ID/msa_default/syn/output"
bash_script_home <- "analysis/code/ID/msa_default"
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
