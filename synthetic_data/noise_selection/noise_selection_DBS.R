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
real_exposure_dbs_all <- PCAWG7::exposure$PCAWG$DBS78
cancer_types <- c(
  "Breast-AdenoCA", "ColoRect-AdenoCA", "Eso-AdenoCA",
  "Kidney-RCC", "Liver-HCC", "Lung-AdenoCA",
  "Ovary-AdenoCA", "Skin-Melanoma", "Stomach-AdenoCA"
)

retval <- lapply(cancer_types, FUN = function(cancer_type) {
  sample_names <-
    grep(
      pattern = cancer_type,
      x = colnames(real_exposure_dbs_all), value = TRUE
    )
  return(real_exposure_dbs_all[, sample_names, drop = FALSE])
})
real_exposure_dbs <- do.call("cbind", retval)

sigs_dbs <- cosmicsig::COSMIC_v3.4$signature$GRCh37$DBS78

# Get the real tumor spectra from selected cancer types
real_spectra_dbs <-
  PCAWG7::spectra$PCAWG$DBS78[, colnames(real_exposure_dbs)]

real_distance_dbs <- get_distance(
  spectra = real_spectra_dbs,
  exposure = real_exposure_dbs,
  sigs = sigs_dbs,
  group = "real"
)

# Add noise to noiseless synthetic data with different negative-binomial size
# parameter
dir_noiseless_dbs <-
  "synthetic_data/DBS/intermed_results/DBS78.syn.exposures.no.noise/"
noiseless_data_dbs <- get_syn_data_info(dir = dir_noiseless_dbs)
seed <- 658220

nb_sizes_dbs <- 1:100

total_cores <- parallel::detectCores()
cores_to_use <- round(total_cores / 2)

noise_data_dbs <- generate_noisy_data(
  seed = seed,
  exposure = noiseless_data_dbs$exposure,
  sigs = sigs_dbs,
  nb_sizes = nb_sizes_dbs,
  mc_cores = cores_to_use
)
syn_distance_dbs <-
  get_multiple_syn_distances(
    list_of_syn_data = noise_data_dbs,
    mc_cores = cores_to_use
  )

p_values <- list()

for (group in names(syn_distance_dbs)) {
  real_distances <- real_distance_dbs[["cosine"]]
  syn_distances <- syn_distance_dbs[[group]][["cosine"]]
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
    "DBS_noise_level_p_value.csv"
  )
)

all_distance_dbs <-
  do.call(dplyr::bind_rows, c(
    list(real_distance_dbs),
    syn_distance_dbs[rownames(df_sorted)[1:5]]
  ))

plot_objects_dbs <-
  create_violinplots(
    distance_df = all_distance_dbs,
    data_type = "DBS", ylim = c(0, 1),
    p_label = "p.format"
  )

ggplot_to_pdf(
  plot_objects = plot_objects_dbs,
  file = file.path(output_home, "noise_selection_DBS.pdf"),
  nrow = 2, ncol = 1,
  width = 8.2677, height = 11.6929, units = "in"
)
