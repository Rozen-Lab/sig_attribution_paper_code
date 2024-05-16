source("analysis/code/analysis_utils.R")

dataset_name <- "SBS"
input_dir <- "analysis/raw_output/SBS/msa_opt/syn/input"
data_identifier <- "syn_data"
output_home <- "analysis/raw_output/SBS/msa_opt/syn/output"
bash_script_home <- "analysis/code/SBS/msa_opt/"

# In this nextflow script, the default values for params.weak_thresholds has been
# changed to ['0.0000', '0.0100', '0.0200', '0.0300']
# See
# https://gitlab.com/s.senkin/MSA/-/wikis/Optimisation-tips#penalty-ranges-to-scan
# for more details
msa_nextflow_file <- 
  "common_code/MSA/run_auto_optimised_analysis_SBS.nf"
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
