if(!grepl("attribution_paper_code", basename(getwd()))) {
  stop("Run this script from the top level directory")
}
rm(list = ls())
message("ls() = ", ls())
source("analysis/code/powerset_analysis.R")

enumeration_of_attributions <- # Used for PCAWG stomach tumors
  function(spectra, 
           selected_sig_names,
           cosine_threshold
  ) {
    
    retval =
      mclapply(
        1:ncol(spectra),
        function(x) {
          cat(x, " ")
          explore_powerset_and_print(
            test_spectrum = spectra[, x, drop = FALSE],
            all_sig_names = selected_sig_names,
            cosine_threshold = cosine_threshold,
            xlabels = FALSE,
            dir_name = ".", # Not used if !plot_pdf
            plot_pdf = FALSE # Too many big files
          )
        },
        # mc.cores = 1) 
        mc.cores = 75) # mclapply
    names(retval) = colnames(spectra)
    return(retval)
  }


# These are the signatures in stomach adenocarcinoma
all_stomach <- mSigAct::ExposureProportions(
  mutation.type = "SBS96",
  cancer.type = "Stomach-AdenoCA"
)
selected_sig_names = names(all_stomach)
# These are 20 signatures, all of which appear in the PCAWG
# stomach spectra

all_stomach_spectra <- pcawg_stomach_spectra()

# If there is no github large file system:
# bigfile = "~/bigdata/powerset_real_spectra.rdata" # PCAWG stomach spectra

# If there *is* github large files system:
bigfile = "output_for_paper/powerset_real_spectra.rdata" # PCAWG stomach spectra
message("Starting to analyze all spectra")
powerset_real_spectra <-
  enumeration_of_attributions(
    spectra = all_stomach_spectra, # [ , 1:10], # for testing
    selected_sig_names = selected_sig_names,

    cosine_threshold = 0.969,
    
  )
save(powerset_real_spectra, file = bigfile)
