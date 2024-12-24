# This file is sourced from the files analysis/code/<mut_type>/run_all_<mut_type>.R 
# stopifnot(basename(getwd()) == "sig_attribution_paper_code")
stopifnot(Sys.getenv("mut_type") != "")
stopifnot(Sys.getenv("sigpro_context_type") != "")
warning("This script does not run MSA")

######

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
for (rare_sig_threshold in global_fitms_rare_sig_thresh) {
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
sigpro_args$seed_in_use = global_random_seed
sigpro_args$python_bin = path.expand("~/software/miniconda3/bin/python")
# sigpro_args$python_bin = "~/miniconda3/envs/sigpro/bin/python"
run_sigpro(Sys.getenv("mut_type"), sigpro_args)

#######

rm(list = ls())
source("analysis/code/yapsa_analysis.R")
# in_per_sample_cutoff of 0.06 based on Vignette, 
# https://bioconductor.org/packages/release/bioc/vignettes/YAPSA/inst/doc/YAPSA.html
# 3.3.4
run_yapsa("SBS", 0.06)

#######

rm(list = ls())
source("analysis/code/deconstruct_analysis.R")
# 0.06 is the default for signature.cutoff
run_deconstruct(Sys.getenv("mut_type"), 0.06)

#######

rm(list = ls())
source("analysis/code/sigest_analysis.R")
run_sigest(Sys.getenv("mut_type"))

#######
if (FALSE) {
  rm(list = ls())
  source("analysis/code/musical_analysis.R")
  run_musical(Sys.getenv("mut_type"), list())
}

#######

rm(list = ls())
source("analysis/code/mutsig_analysis.R")
run_mutsig(Sys.getenv("mut_type"))

#######

if (Sys.getenv("mut_type") == "SBS") {
  rm(list = ls())
  source("analysis/code/siglasso_analysis.R")
  run_siglasso(Sys.getenv("mut_type"), list(use_prior = FALSE))
}

#######

if (Sys.getenv("mut_type") == "SBS") {
  rm(list = ls())
  source("analysis/code/siglasso_analysis.R")
  run_siglasso(Sys.getenv("mut_type"), more_args = list(use_prior = TRUE))
}

#######

rm(list = ls())
source("analysis/code/sigfit_analysis.R")
run_sigfit("SBS")

#######

rm(list = ls())
source("analysis/code/sigspack_analysis.R")
run_sigspack("SBS")

