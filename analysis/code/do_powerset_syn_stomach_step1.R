if(!grepl("attribution_paper_code", basename(getwd()))) {
  stop("Run this script from the top level directory")
}
rm(list = ls())
library(writexl)
source("analysis/code/powerset_analysis.R")
source("analysis/code/reoptimize_similarities_from_exposures.R")


get_alternative_attribution <-
  function(sample_names, spectra, selected_sig_names) {
    # browser()
    mc.cores = 100
    if (Sys.info()["sysname"] == "Windows") mc_cores = 1
    stop("Relace this with explore_powerset_noprint")
    
    retval <- parallel::mclapply(
      sample_names,
      function(x) {
        cat(x, "\n")
        ret <-
          list(
            sample = x,
            explore_powerset_noprint(
              test_spectrum = spectra[, x, drop = FALSE],
              all_sig_names = selected_sig_names,
              cosine_threshold = 0
            )
          )
        return(ret)
      },
      mc.cores = mc.cores
    )

    return(retval)
  }

syn_stomach_spectra = synthetic_stomach_spectra()
syn_stomach_exposures = synthetic_stomach_exposures()
sigs_to_use = rownames(syn_stomach_exposures)[rowSums(syn_stomach_exposures) > 0]

sample_names = colnames(syn_stomach_spectra) # for testing

syn_stomach_powerset <-
  get_alternative_attribution(
    sample_names = sample_names,
    spectra = syn_stomach_spectra,
    selected_sig_names = sigs_to_use
  )

save(syn_stomach_powerset, 
     # If no git LFS
     # file = "~/bigdata/syn_stomach_powerset.rdata"
     # If yes git LFS
     file = "output_for_paper/syn_stomach_powerset.rdata"
)

