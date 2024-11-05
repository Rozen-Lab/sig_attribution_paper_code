rm(list = ls())

source("synthetic_data/data_gen_code/data_gen_utils.R")

library(dplyr)
library(ggpubr)
library(gridExtra)
library(SynSigGen) # remotes::install_github(repo = "steverozen/SynSigGen",ref = "1.2.2-branch")
library(ICAMS) # remotes::install_github("steverozen/ICAMS", ref = "v3.0.8-branch")
library(mSigAct) # remotes::install_github(repo = "steverozen/mSigAct", ref = "v3.0.1-branch")
library(PCAWG7) # remotes::install_github(repo = "steverozen/PCAWG7", ref = "v0.1.3-branch")
library(cosmicsig) # remotes::install_github(repo = "Rozen-Lab/cosmicsig", ref = "v1.2.0-branch")

# Get the real exposure for tumors from selected cancer types
real_exposure_sbs_all <- PCAWG7::exposure$PCAWG$SBS96
cancer_types <- c(
  "Breast-AdenoCA", "ColoRect-AdenoCA", "Eso-AdenoCA",
  "Kidney-RCC", "Liver-HCC", "Lung-AdenoCA",
  "Ovary-AdenoCA", "Skin-Melanoma", "Stomach-AdenoCA"
)

retval <- lapply(cancer_types, FUN = function(cancer_type) {
  sample_names <-
    grep(
      pattern = cancer_type,
      x = colnames(real_exposure_sbs_all), value = TRUE
    )
  return(real_exposure_sbs_all[, sample_names, drop = FALSE])
})
real_exposure_sbs <- do.call("cbind", retval)
rownames(real_exposure_sbs) <-
  gsub(
    pattern = "SBS22", replacement = "SBS22a",
    x = rownames(real_exposure_sbs)
  )

sigs_sbs <- cosmicsig::COSMIC_v3.4$signature$GRCh37$SBS96
sigs_sbs <-
  sigs_sbs[, setdiff(
    colnames(sigs_sbs),
    paste0("SBS40", letters[1:3])
  )]
sbs40 <- cosmicsig::COSMIC_v3.3$signature$GRCh37$SBS96[, "SBS40", drop = FALSE]
sigs_sbs_v2 <- cbind(sigs_sbs, sbs40)

# Get the real tumor spectra from selected cancer types
real_spectra_sbs <-
  PCAWG7::spectra$PCAWG$SBS96[, colnames(real_exposure_sbs)]

real_distance_sbs <- get_distance(
  spectra = real_spectra_sbs,
  exposure = real_exposure_sbs,
  sigs = sigs_sbs_v2,
  group = "real"
)

# Add noise to noiseless synthetic data with different negative-binomial size
# parameter
dir_noiseless_sbs <-
  "synthetic_data/SBS/intermed_results/SBS96.syn.exposures.no.noise/"
noiseless_data_sbs <- get_syn_data_info(dir = dir_noiseless_sbs)
seed <- 658220

nb_sizes_sbs <- 1:100

total_cores <- parallel::detectCores()
cores_to_use <- round(total_cores / 2)

sigs_sbs_v3 <- 
  cbind(cosmicsig::COSMIC_v3.4$signature$GRCh37$SBS96, sbs40)

noise_data_sbs <- generate_noisy_data(
  seed = seed,
  exposure = noiseless_data_sbs$exposure,
  sigs = sigs_sbs_v3,
  nb_sizes = nb_sizes_sbs,
  mc_cores = cores_to_use
)
syn_distance_sbs <-
  get_multiple_syn_distances(
    list_of_syn_data = noise_data_sbs,
    mc_cores = cores_to_use
  )

p_values <- list()

for (group in names(syn_distance_sbs)) {
  real_distances <- real_distance_sbs[["cosine"]]
  syn_distances <- syn_distance_sbs[[group]][["cosine"]]
  test_results <- stats::wilcox.test(x = real_distances, y = syn_distances)
  p_values[[group]] <- test_results$p.value
}

output_home <- "synthetic_data/noise_selection/output"

if (!dir.exists(output_home)) {
  dir.create(output_home, recursive = TRUE)
}

df <- data.frame(p_value = unlist(p_values))
df_sorted <- dplyr::arrange(df, dplyr::desc(p_value))
write.csv(df_sorted,
  file = file.path(
    output_home,
    "SBS_noise_level_p_value.csv"
  )
)

all_distance_sbs <-
  do.call(dplyr::bind_rows, c(
    list(real_distance_sbs),
    syn_distance_sbs[rownames(df_sorted)[1:5]]
  ))

plot_objects_sbs <-
  create_violinplots(
    distance_df = all_distance_sbs,
    data_type = "SBS", ylim = c(0, 1),
    p_label = "p.format"
  )

ggplot_to_pdf(
  plot_objects = plot_objects_sbs,
  file = file.path(output_home, "noise_selection_SBS.pdf"),
  nrow = 2, ncol = 1,
  width = 8.2677, height = 11.7929, units = "in"
)
