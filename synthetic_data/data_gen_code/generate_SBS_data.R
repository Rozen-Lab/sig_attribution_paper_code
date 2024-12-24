# Run this script from the top level directory

source("synthetic_data/data_gen_code/data_gen_utils.R")

library(cosmicsig)
library(PCAWG7)
library(ICAMS)
library(SynSigGen)
library(mSigTools)

##################################################################
##                      Data preprocessing                      ##
##################################################################
# Get the real exposures from PCAWG assignments
real_exposures_sbs96 <- PCAWG7::exposure$PCAWG$SBS96
pcawg_sbs96_catalog <- PCAWG7::spectra$PCAWG$SBS96
rownames(real_exposures_sbs96) <-
  gsub(
    pattern = "SBS22", replacement = "SBS22a",
    x = rownames(real_exposures_sbs96)
  )

# Only select samples that belong to the selected cancer types
cancer_types <- c(
  "Breast-AdenoCA", "ColoRect-AdenoCA", "Eso-AdenoCA",
  "Liver-HCC", "Lung-AdenoCA",
  "Ovary-AdenoCA", "Skin-Melanoma", "Stomach-AdenoCA"
)
indices_selected_types <- unlist(sapply(cancer_types, FUN = function(x) {
  grep(x, colnames(real_exposures_sbs96))
}))
real_exposures_sbs96 <- real_exposures_sbs96[, indices_selected_types]

# Exclude samples which have mutations less than 100
samples_less_than_100 <- names(which(colSums(pcawg_sbs96_catalog) < 100))
indices_less_than_100 <-
  which(colnames(real_exposures_sbs96) %in% samples_less_than_100)
real_exposures_sbs96 <-
  real_exposures_sbs96[, -indices_less_than_100, drop = FALSE]

# Exclude tumors which have possible artifact signatures
artifact_sigs <- cosmicsig::possible_artifacts()
indices_artifact_sigs <-
  which(rownames(real_exposures_sbs96) %in% artifact_sigs)
artifact_sigs_selected_types <-
  rownames(real_exposures_sbs96)[indices_artifact_sigs]
tumors_to_remove <- sapply(artifact_sigs_selected_types, FUN = function(x) {
  exposure <- real_exposures_sbs96[x, ]
  return(names(exposure[exposure > 0]))
})
tumors_to_remove <- unique(unlist(tumors_to_remove))
real_exposures_sbs96 <-
  real_exposures_sbs96[, !colnames(real_exposures_sbs96) %in% tumors_to_remove]

real_exposures_sbs96 <- remove_zero_activity_sigs(real_exposures_sbs96)

# Get the real kidney exposures from Senkin et al., 2024
# https://doi.org/10.1038/s41586-024-07368-2
kidney_exp_file <-
  "synthetic_data/data_gen_code/Senkin_et_al.,_2024_SBS_exposures.csv"
real_kidney_exposures_sbs96 <-
  mSigTools::read_exposure(file = kidney_exp_file)
colnames(real_kidney_exposures_sbs96) <-
  paste0("Kidney-RCC::", colnames(real_kidney_exposures_sbs96))

sum(colSums(real_kidney_exposures_sbs96) < 100)
sum(rownames(real_kidney_exposures_sbs96) %in% artifact_sigs)

##################################################################
##   Calculate number of synthetic tumors in each cancer type   ##
##################################################################

msi_sample_indices_selected_types <-
  unlist(sapply(pcawg_msi_tumor_ids, FUN = function(x) {
    grep(x, colnames(real_exposures_sbs96))
  }))
msi_sample_ids <- names(msi_sample_indices_selected_types)
length(msi_sample_ids) # 19

pole_sample_indices_selected_types <-
  unlist(sapply(pcawg_pole_tumor_ids, FUN = function(x) {
    grep(x, colnames(real_exposures_sbs96))
  }))
pole_sample_ids <- names(pole_sample_indices_selected_types)
length(pole_sample_ids) # 8

real_exposures_sbs96_no_msi_pole <-
  real_exposures_sbs96[, -c(
    msi_sample_indices_selected_types,
    pole_sample_indices_selected_types
  ), drop = FALSE]
real_exposures_sbs96_no_msi_pole <-
  remove_zero_activity_sigs(real_exposures_sbs96_no_msi_pole)
real_exposures_sbs96_msi <-
  real_exposures_sbs96[, msi_sample_indices_selected_types, drop = FALSE]
real_exposures_sbs96_msi <- remove_zero_activity_sigs(real_exposures_sbs96_msi)
real_exposures_sbs96_pole <-
  real_exposures_sbs96[, pole_sample_indices_selected_types, drop = FALSE]
real_exposures_sbs96_pole <- remove_zero_activity_sigs(real_exposures_sbs96_pole)

num_samples_total <- calculate_num_samples(real_exposures_sbs96)
num_samples_msi <- calculate_num_samples(real_exposures_sbs96_msi)
cancer_types_msi <- names(num_samples_msi)
num_samples_pole <- calculate_num_samples(real_exposures_sbs96_pole)
cancer_types_pole <- names(num_samples_pole)
num_samples_no_msi_pole <- calculate_num_samples(real_exposures_sbs96_no_msi_pole)

# Generate 100 synthetic tumors for the selected cancer types (total 800). Scale
# the original number of tumors in each cancer type in real exposure accordingly
scale_factors <- 100 / num_samples_total

# Calculate the number of MSI-H synthetic tumors in each cancer type
num_samples_msi_scaled <-
  sapply(cancer_types_msi, FUN = function(x) {
    scaled_number <- ceiling(scale_factors[x] * num_samples_msi[x])
    names(scaled_number) <- NULL
    return(scaled_number)
  })

# Calculate the number of POLE-mutated synthetic tumors in each cancer type
num_samples_pole_scaled <-
  sapply(cancer_types_pole, FUN = function(x) {
    scaled_number <- ceiling(scale_factors[x] * num_samples_pole[x])
    names(scaled_number) <- NULL
    return(scaled_number)
  })

# Calculate the number of non MSI-H, non POLE-mutated synthetic tumors in each
# cancer type
num_samples_no_msi_pole_scaled <- rep(100, length(cancer_types))
names(num_samples_no_msi_pole_scaled) <- cancer_types
for (i in cancer_types_msi) {
  num_samples_no_msi_pole_scaled[i] <-
    num_samples_no_msi_pole_scaled[i] - num_samples_msi_scaled[i]
}
num_samples_no_msi_pole_scaled[cancer_types_pole] <-
  num_samples_no_msi_pole_scaled[cancer_types_pole] - num_samples_pole_scaled

# Make sure the total number of synthetic tumors is 800
sum(num_samples_msi_scaled) + sum(num_samples_pole_scaled) +
  sum(num_samples_no_msi_pole_scaled)

# Calculate number of MSI-H samples in real kidney
msi_sig_names <-
  c("SBS6", "SBS14", "SBS15", "SBS20", "SBS21", "SBS26", "SBS44")
tmp <-
  real_kidney_exposures_sbs96[rownames(real_kidney_exposures_sbs96) %in%
    msi_sig_names, ]
tmp2 <- tmp[, colSums(tmp) > 0]

kidney_msi_sample_ids <- colnames(tmp2)
length(kidney_msi_sample_ids) # 6

kidney_msi_sample_indices <-
  unlist(sapply(kidney_msi_sample_ids, FUN = function(x) {
    grep(x, colnames(real_kidney_exposures_sbs96))
  }))
real_kidney_msi_exposures_sbs96 <-
  real_kidney_exposures_sbs96[, kidney_msi_sample_indices, drop = FALSE]
real_kidney_msi_exposures_sbs96 <-
  remove_zero_activity_sigs(real_kidney_msi_exposures_sbs96)

real_kidney_exposures_sbs96_no_msi <-
  real_kidney_exposures_sbs96[, -kidney_msi_sample_indices, drop = FALSE]
real_kidney_exposures_sbs96_no_msi <-
  remove_zero_activity_sigs(real_kidney_exposures_sbs96_no_msi)

num_kidney_samples_total <- ncol(real_kidney_exposures_sbs96)
num_kidney_samples_msi <- ncol(real_kidney_msi_exposures_sbs96)
num_kidney_samples_no_msi <- ncol(real_kidney_exposures_sbs96_no_msi)

# Only generate 100 synthetic tumors for kidney. Scale the original number of
# tumors in real exposure accordingly
kidney_scale_factor <- 100 / num_kidney_samples_total

# Calculate the number of MSI-H (MSI-High) synthetic tumors in kidney
num_samples_kidney_msi_scaled <-
  ceiling(kidney_scale_factor * num_kidney_samples_msi)

# Calculate the number of non MSI-H synthetic tumors in kidney
num_samples_kidney_no_msi_scaled <- 100 - num_samples_kidney_msi_scaled

##################################################################
##                 Generation of synthetic data                 ##
##################################################################

output_dir_sbs96_no_msi_pole <-
  "./synthetic_data/SBS/intermed_results/PCAWG.SBS96.syn.exposures.no.msi"
output_dir_sbs96_msi <-
  "./synthetic_data/SBS/intermed_results/PCAWG.SBS96.syn.exposures.msi"
output_dir_sbs96_pole <-
  "./synthetic_data/SBS/intermed_results/PCAWG.SBS96.syn.exposures.pole"

distribution <- "neg.binom"
sample_prefix_name <- ""
mutation_type <- "SBS96"
seed <- 658220
input_sigs_sbs96 <- cosmicsig::COSMIC_v3.4$signature$GRCh37$SBS96
input_sigs_sbs96 <-
  input_sigs_sbs96[, setdiff(
    colnames(input_sigs_sbs96),
    paste0("SBS40", letters[1:3])
  )]
sbs40 <- cosmicsig::COSMIC_v3.3$signature$GRCh37$SBS96[, "SBS40", drop = FALSE]
input_sigs_sbs96 <- cbind(input_sigs_sbs96, sbs40)

real_exposures_sbs96_all <-
  SynSigGen::MergeExposures(
    list.of.exposures =
      list(
        real_exposures_sbs96,
        real_kidney_exposures_sbs96
      )
  )

sig_params_sbs96_all_types <-
  SynSigGen:::GetSynSigParamsFromExposures(
    exposures = real_exposures_sbs96_all,
    distribution = distribution,
    sig.params = SynSigGen::signature.params$SBS96
  )

# Generate synthetic tumors that are not MSI-H or POLE-mutated
synthetic_tumors_sbs96_no_msi_pole <-
  SynSigGen::GenerateSyntheticTumors(
    seed = seed,
    dir = output_dir_sbs96_no_msi_pole,
    cancer.types = cancer_types,
    samples.per.cancer.type = num_samples_no_msi_pole_scaled,
    input.sigs = input_sigs_sbs96,
    real.exposures = real_exposures_sbs96_no_msi_pole,
    distribution = distribution,
    sample.prefix.name = sample_prefix_name,
    sig.params = sig_params_sbs96_all_types
  )
unlink(output_dir_sbs96_no_msi_pole, recursive = TRUE)
syn_exposures_sbs96_no_msi_pole <-
  synthetic_tumors_sbs96_no_msi_pole$ground.truth.exposures
syn_exposures_sbs96_no_msi_pole <-
  check_syn_exposures(syn_exposures_sbs96_no_msi_pole)

# Generate MSI-H synthetic tumors
synthetic_tumors_sbs96_msi <-
  generate_subtype_syn_tumors(
    seed = seed,
    dir = output_dir_sbs96_msi,
    cancer_types = cancer_types_msi,
    samples_per_caner_type = num_samples_msi_scaled,
    input_sigs = input_sigs_sbs96,
    real_exposure = real_exposures_sbs96_msi,
    distribution = distribution,
    sample_prefix_name = sample_prefix_name,
    tumor_marker_name = "MSI-H",
    sig_params = sig_params_sbs96_all_types
  )
unlink(output_dir_sbs96_msi, recursive = TRUE)
syn_exposures_sbs96_msi <-
  synthetic_tumors_sbs96_msi$ground.truth.exposures

# Generate POLE-mutated synthetic tumors
synthetic_tumors_sbs96_pole <-
  generate_subtype_syn_tumors(
    seed = seed,
    dir = output_dir_sbs96_pole,
    cancer_types = cancer_types_pole,
    samples_per_caner_type = num_samples_pole_scaled,
    input_sigs = input_sigs_sbs96,
    real_exposure = real_exposures_sbs96_pole,
    distribution = distribution,
    sample_prefix_name = sample_prefix_name,
    tumor_marker_name = "POLE",
    sig_params = sig_params_sbs96_all_types
  )
unlink(output_dir_sbs96_pole, recursive = TRUE)
syn_exposures_sbs96_pole <-
  synthetic_tumors_sbs96_pole$ground.truth.exposures

# Combine the non MSI-H non POLE, MSI-H and POLE synthetic exposures in each
# cancer type
synthetic_exposures_sbs96 <-
  combine_exposure(
    syn_exposures_sbs96_no_msi_pole,
    syn_exposures_sbs96_pole,
    syn_exposures_sbs96_msi
  )

# Generation of synthetic kidney data
output_dir_kidney_sbs96_no_msi <-
  "./synthetic_data/SBS/intermed_results/kidney.SBS96.syn.exposures.no.msi"
output_dir_kidney_sbs96_msi <-
  "./synthetic_data/SBS/intermed_results/kidney.SBS96.syn.exposures.msi"

input_sigs_sbs96_v2 <- cosmicsig::COSMIC_v3.4$signature$GRCh37$SBS96

# Generate non MSI-H synthetic kidney tumors
synthetic_tumors_kidney_sbs96_no_msi <-
  SynSigGen::GenerateSyntheticTumors(
    seed = seed,
    dir = output_dir_kidney_sbs96_no_msi,
    cancer.types = "Kidney-RCC",
    samples.per.cancer.type = num_samples_kidney_no_msi_scaled,
    input.sigs = input_sigs_sbs96_v2,
    real.exposures = real_kidney_exposures_sbs96_no_msi,
    distribution = distribution,
    sample.prefix.name = sample_prefix_name,
    sig.params = sig_params_sbs96_all_types
  )
unlink(output_dir_kidney_sbs96_no_msi, recursive = TRUE)
syn_exposures_kidney_sbs96_no_msi <-
  synthetic_tumors_kidney_sbs96_no_msi$ground.truth.exposures

# Generate MSI-H synthetic kidney tumors
synthetic_tumors_kidney_sbs96_msi <-
  generate_subtype_syn_tumors(
    seed = seed,
    dir = output_dir_kidney_sbs96_msi,
    cancer_types = "Kidney-RCC",
    samples_per_caner_type = num_samples_kidney_msi_scaled,
    input_sigs = input_sigs_sbs96_v2,
    real_exposure = real_kidney_msi_exposures_sbs96,
    distribution = distribution,
    sample_prefix_name = sample_prefix_name,
    tumor_marker_name = "MSI-H",
    sig_params = sig_params_sbs96_all_types
  )
unlink(output_dir_kidney_sbs96_msi, recursive = TRUE)
syn_exposures_kidney_sbs96_msi <-
  synthetic_tumors_kidney_sbs96_msi$ground.truth.exposures

# Combine the non MSI-H and MSI-H synthetic exposures in kidney
synthetic_exposures_kidney_sbs96 <-
  combine_exposure(
    syn_exposures_kidney_sbs96_no_msi,
    syn_exposures_kidney_sbs96_msi
  )

# Generate the combined synthetic tumors
synthetic_exposures_sbs96_all <-
  SynSigGen::MergeExposures(
    list.of.exposures =
      list(
        synthetic_exposures_sbs96,
        synthetic_exposures_kidney_sbs96
      )
  )

output_dir_sbs96 <-
  "./synthetic_data/SBS/intermed_results/SBS96.syn.exposures.no.noise"
output_dir_sbs96_nb_size_selected <-
  "./synthetic_data/SBS/intermed_results/SBS96.syn.exposures.noisy.neg.binom.size.selected"

write_sig_params(
  dir = output_dir_sbs96,
  real_exposure = real_exposures_sbs96_all,
  synthetic_exposure = synthetic_exposures_sbs96_all,
  cancer_types = c(cancer_types, "Kidney-RCC"),
  distribution = distribution,
  sig_params = sig_params_sbs96_all_types,
  sample_prefix_name = sample_prefix_name,
  mutation_type = mutation_type
)

input_sigs_sbs96_v3 <- cbind(input_sigs_sbs96_v2, sbs40)

catalog.info <- SynSigGen:::CreateSynCatalogs(
  signatures = input_sigs_sbs96_v3,
  exposures  = synthetic_exposures_sbs96_all
)
ICAMS::WriteCatalog(
  catalog.info$ground.truth.catalog,
  file = file.path(output_dir_sbs96, "ground.truth.syn.catalog.csv")
)
mSigTools::write_exposure(
  exposure = synthetic_exposures_sbs96_all,
  file     = file.path(output_dir_sbs96, "ground.truth.syn.exposures.csv")
)

# Add noise to the synthetic tumors
sbs96_noisy_tumors_size_selected <-
  SynSigGen::GenerateNoisyTumors(
    seed = seed,
    dir = output_dir_sbs96_nb_size_selected,
    input.exposure = synthetic_exposures_sbs96_all,
    signatures = input_sigs_sbs96_v3,
    n.binom.size = 9
  )

noisy_exposures_size_selected_sbs96 <- 
  sbs96_noisy_tumors_size_selected$exposures

get_sig_universes(
  exposures = real_exposures_sbs96_all,
  filename = file.path(
    output_dir_sbs96_nb_size_selected,
    "ground.truth.sig.universe.csv"
  ), # This is signature names for each cancer type
  sigs = input_sigs_sbs96_v3,
  sig_file = file.path(
    output_dir_sbs96_nb_size_selected,
    "ground.truth.sigs.csv"
  )
)

#################################################################
##                   Plot data distributions                   ##
#################################################################

data_distribution_file <-
  "synthetic_data/SBS/intermed_results/SBS_syn_data_distribution.pdf"
grDevices::pdf(
  file = data_distribution_file,
  width = 11.6929, height = 8.2677, onefile = TRUE
)
par(mfrow = c(3, 3))
plot_exposure_distribution(
  real_exposure = real_exposures_sbs96_all,
  synthetic_exposure = synthetic_exposures_sbs96_all,
  noisy_exposure = noisy_exposures_size_selected_sbs96,
  size = 9,
  distribution = distribution,
  sig_params = sig_params_sbs96_all_types,
  sample_prefix_name = sample_prefix_name
)
grDevices::dev.off()


#################################################################
##               Put data in standard places / file names      ##
#################################################################

source("synthetic_data/data_gen_code/data_gen_rename.R")

old_dataset_names <- basename(output_dir_sbs96_nb_size_selected)

data_gen_rename(
  dataset = "SBS",
  old_dataset_name = old_dataset_names,
  top_folder_name = "synthetic_data"
)
