# Run this script with the top level directory as the working directory
# 



stopifnot(Sys.getenv("mut_type") != "")
stopifnot(Sys.getenv("sigpro_context_type") != "")


######

stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/pasa_analysis.R")

total_cores <- parallel::detectCores() # 256 on the machine we are using
mc_cores_per_sample <- min(5, total_cores)
pasa_args = list()
pasa_args$mc_cores_per_sample  <- mc_cores_per_sample
pasa_args$num_parallel_samples <- floor(total_cores / mc_cores_per_sample)
pasa_args$seed_in_use = seed = 145879
rm(mc_cores_per_sample, total_cores)
run_pasa(mut_type = Sys.getenv("mut_type"),  more_args = pasa_args)
rm(pasa_args)

#######

rm(list = ls())
source("analysis/code/fitms_analysis.R")
for (rare_sig_threshold in c(0.01, 0.03, 0.06)) {
  run_fitms(Sys.getenv("mut_type"), rare_sig_threshold)
}

#######

rm(list = ls())
source("analysis/code/mp_analysis.R")
run_mp(Sys.getenv("mut_type"))

#######

rm(list = ls())
source("analysis/code/sigpro_analysis.R")
sigpro_args = list()
sigpro_args$context_type = Sys.getenv("sigpro_context_type")
sigpro_args$seed_in_use = 145879
sigpro_args$python_bin = path.expand("~/miniconda3/envs/sigpro/bin/python")
# sigpro_args$python_bin = "/home/e0012078/software/miniconda3/bin/python"
run_sigpro(Sys.getenv("mut_type"), sigpro_args)

#######

rm(list = ls())
source("analysis/code/yapsa_analysis.R")
# default in_per_sample_cutoff is 0
for (in_per_sample_cutoff in c(0, 0.01, 0.03, 0.06, 0.1)) {
  run_yapsa(Sys.getenv("mut_type"), in_per_sample_cutoff)
}

#######

rm(list = ls())
source("analysis/code/yapsa_analysis.R")
# default in_per_sample_cutoff is 0
for (in_per_sample_cutoff in c(0, 0.01, 0.03, 0.06, 0.1)) {
  run_yapsa(Sys.getenv("mut_type"), in_per_sample_cutoff)
}

#######

rm(list = ls())
source("analysis/code/deconstruct_analysis.R")
# 0.06 is the default for signature.cutoff
for (signature.cutoff in c(0, 0.03, 0.06, 0.1)) {
  run_deconstruct(Sys.getenv("mut_type"), signature.cutoff)
}

#######

rm(list = ls())
source("analysis/code/sigest_analysis.R")
run_sigest(Sys.getenv("mut_type"))

