if(basename(getwd()) != "sig_attribution_paper_code") {
  stop("Run this script from the top level directory")
}

source("analysis/code/analysis_utils.R")

# These are a subset of the signatures in stomach adenocarcinoma
all_stomach <- mSigAct::ExposureProportions(
  mutation.type = "SBS96",
  cancer.type = "Stomach-AdenoCA"
)

# Working on the power set of all signatures in all_stomach is too slow.
# Get rid of the signatures seen less frequently.
selected_sig_names <- names(all_stomach[all_stomach > min(all_stomach)])

s96 <- PCAWG7::spectra$PCAWG$SBS96
all_stomach_spectra <-
  s96[, grep("Stomach-AdenoCA", colnames(s96), fixed = TRUE)]

# Exhaustive enumeration of stomach cancer attributions
if (!exists("explore_real_spectra")) {
  if (file.exists("common_code/figure_2/explore_real_spectra.Rdata")) {
    load("common_code/figure_2/explore_real_spectra.Rdata")
  } else {
    dir_name <- "common_code/figure_2/enumeration_of_attributions"
    explore_real_spectra <-
      enumeration_of_attributions(
        spectra = all_stomach_spectra,
        selected_sig_names = selected_sig_names,
        dir_name = dir_name
      )
    save(explore_real_spectra,
      file = "common_code/figure_2/explore_real_spectra.Rdata"
    )
  }
}

all_exp <- lapply(explore_real_spectra, `[[`, "exposures")
num_recon <- unlist(lapply(all_exp, number_reconstructions))

sum(num_recon > 0)
# 49
median(num_recon[num_recon > 0])
# 257

file <- "common_code/figure_2/figure_2_panel_A.pdf"
figure_2_histogram(
  file = file,
  num_of_recon = num_recon[num_recon > 0]
)

# Generate main text figure showing manually selected example reconstructions
# and generate corresponding draft supplementary table
sample_name <- "Stomach-AdenoCA::SP85251"
figure_2_and_sup_table_s2(
  sample_name = sample_name,
  selected_sig_names = selected_sig_names,
  dir_name = "common_code/figure_2/"
)
