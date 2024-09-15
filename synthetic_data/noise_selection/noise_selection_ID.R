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
real_exposure_id_all <- PCAWG7::exposure$PCAWG$ID
cancer_types <- c(
  "Breast-AdenoCA", "ColoRect-AdenoCA", "Eso-AdenoCA",
  "Kidney-RCC", "Liver-HCC", "Lung-AdenoCA",
  "Ovary-AdenoCA", "Skin-Melanoma", "Stomach-AdenoCA"
)

retval <- lapply(cancer_types, FUN = function(cancer_type) {
  sample_names <-
    grep(
      pattern = cancer_type,
      x = colnames(real_exposure_id_all), value = TRUE
    )
  return(real_exposure_id_all[, sample_names, drop = FALSE])
})
real_exposure_id <- do.call("cbind", retval)

sigs_id <- cosmicsig::COSMIC_v3.4$signature$GRCh37$ID

# Get the real tumor spectra from selected cancer types
real_spectra_id <-
  PCAWG7::spectra$PCAWG$ID[, colnames(real_exposure_id)]

real_distance_id <- get_distance(
  spectra = real_spectra_id,
  exposure = real_exposure_id,
  sigs = sigs_id,
  group = "real"
)

# Add noise to noiseless synthetic data with different negative-binomial size
# parameter
dir_noiseless_id <-
  "synthetic_data/ID/intermed_results/ID.syn.exposures.no.noise/"
noiseless_data_id <- get_syn_data_info(dir = dir_noiseless_id)
seed <- 658220

nb_sizes_id <- 1:100

total_cores <- parallel::detectCores()
cores_to_use <- round(total_cores / 2)

noise_data_id <- generate_noisy_data(
  seed = seed,
  exposure = noiseless_data_id$exposure,
  sigs = sigs_id,
  nb_sizes = nb_sizes_id,
  mc_cores = cores_to_use
)
syn_distance_id <-
  get_multiple_syn_distances(
    list_of_syn_data = noise_data_id,
    mc_cores = cores_to_use
  )

p_values <- list()

for (group in names(syn_distance_id)) {
  real_distances <- real_distance_id[["cosine"]]
  syn_distances <- syn_distance_id[[group]][["cosine"]]
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
    "ID_noise_level_p_value.csv"
  )
)

all_distance_id <-
  do.call(dplyr::bind_rows, c(
    list(real_distance_id),
    syn_distance_id[rownames(df_sorted)[1:5]]
  ))

plot_objects_id <-
  create_violinplots(
    distance_df = all_distance_id,
    data_type = "ID", ylim = c(0, 1),
    p_label = "p.format"
  )

ggplot_to_pdf(
  plot_objects = plot_objects_id,
  file = file.path(output_home, "noise_selection_ID.pdf"),
  nrow = 2, ncol = 1,
  width = 8.2677, height = 11.6929, units = "in"
)
