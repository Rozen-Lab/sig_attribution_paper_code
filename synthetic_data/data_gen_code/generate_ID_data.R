# Run this script from the top level directory

source("synthetic_data/data_gen_code/data_gen_utils.R")

library(cosmicsig)
library(PCAWG7)
library(SynSigGen)

##################################################################
##                      Data preprocessing                      ##
##################################################################
# Get the real exposures from PCAWG assignments
real_exposures_id <- PCAWG7::exposure$PCAWG$ID
pcawg_id_catalog <- PCAWG7::spectra$PCAWG$ID

# Only select samples that belong to the selected cancer types
cancer_types <- c(
  "Breast-AdenoCA", "ColoRect-AdenoCA", "Eso-AdenoCA",
  "Kidney-RCC", "Liver-HCC", "Lung-AdenoCA",
  "Ovary-AdenoCA", "Skin-Melanoma", "Stomach-AdenoCA"
)
indices_nine_types <- unlist(sapply(cancer_types, FUN = function(x) {
  grep(x, colnames(real_exposures_id))
}))
real_exposures_id <- real_exposures_id[, indices_nine_types]

# Exclude samples which have mutations less than 100
samples_less_than_100 <- names(which(colSums(pcawg_id_catalog) < 100))
indices_less_than_100 <-
  which(colnames(real_exposures_id) %in% samples_less_than_100)
real_exposures_id <-
  real_exposures_id[, -indices_less_than_100, drop = FALSE]

real_exposures_id <- remove_zero_activity_sigs(real_exposures_id)

##################################################################
##   Calculate number of synthetic tumors in each cancer type   ##
##################################################################

msi_sample_indices_nine_types <-
  unlist(sapply(pcawg_msi_tumor_ids, FUN = function(x) {
    grep(x, colnames(real_exposures_id))
  }))
msi_sample_ids <- names(msi_sample_indices_nine_types)
length(msi_sample_ids) # 20

pole_sample_indices_nine_types <-
  unlist(sapply(pcawg_pole_tumor_ids, FUN = function(x) {
    grep(x, colnames(real_exposures_id))
  }))
pole_sample_ids <- names(pole_sample_indices_nine_types)

# There are no POLE samples in real_exposures_id
length(pole_sample_ids) # 8

real_exposures_id_no_msi_pole <-
  real_exposures_id[, -c(
    msi_sample_indices_nine_types,
    pole_sample_indices_nine_types
  ), drop = FALSE]
real_exposures_id_no_msi_pole <-
  remove_zero_activity_sigs(real_exposures_id_no_msi_pole)
real_exposures_id_msi <-
  real_exposures_id[, msi_sample_indices_nine_types, drop = FALSE]
real_exposures_id_msi <- remove_zero_activity_sigs(real_exposures_id_msi)
real_exposures_id_pole <-
  real_exposures_id[, pole_sample_indices_nine_types, drop = FALSE]
real_exposures_id_pole <- remove_zero_activity_sigs(real_exposures_id_pole)

num_samples_total <- calculate_num_samples(real_exposures_id)
num_samples_msi <- calculate_num_samples(real_exposures_id_msi)
cancer_types_msi <- names(num_samples_msi)
num_samples_pole <- calculate_num_samples(real_exposures_id_pole)
cancer_types_pole <- names(num_samples_pole)
num_samples_no_msi_pole <- calculate_num_samples(real_exposures_id_no_msi_pole)

# Generate 100 synthetic tumors for the nine cancer types (total 900). Scale
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

# Make sure the total number of synthetic tumors is 270
sum(num_samples_msi_scaled) + sum(num_samples_pole_scaled) +
  sum(num_samples_no_msi_pole_scaled)

##################################################################
##                 Generation of synthetic data                 ##
##################################################################

output_dir_id_no_msi_pole <- "./synthetic_data/ID/intermed_results/ID.syn.exposures.no.msi.pole"
output_dir_id_msi <- "./synthetic_data/ID/intermed_results/ID.syn.exposures.msi"
output_dir_id_pole <- "./synthetic_data/ID/intermed_results/ID.syn.exposures.pole"
output_dir_id <- "./synthetic_data/ID/intermed_results/ID.syn.exposures.no.noise"
output_dir_id_nb_size_selected <-
  "./synthetic_data/ID/intermed_results/ID.syn.exposures.noisy.neg.binom.size.selected"

distribution <- "neg.binom"
sample_prefix_name <- ""
mutation_type <- "ID"
seed <- 658220
input_sigs_id <- cosmicsig::COSMIC_v3.4$signature$GRCh37$ID

sig_params_id_nine_types <-
  SynSigGen:::GetSynSigParamsFromExposures(
    exposures = real_exposures_id,
    distribution = distribution,
    sig.params = SynSigGen::signature.params$ID
  )

# Generate synthetic tumors that are not MSI-H or POLE-mutated
syn_tumors_id_no_msi_pole <-
  SynSigGen::GenerateSyntheticTumors(
    seed = seed,
    dir = output_dir_id_no_msi_pole,
    cancer.types = cancer_types,
    samples.per.cancer.type = num_samples_no_msi_pole_scaled,
    input.sigs = input_sigs_id,
    real.exposures = real_exposures_id_no_msi_pole,
    distribution = distribution,
    sample.prefix.name = sample_prefix_name,
    sig.params = sig_params_id_nine_types
  )
unlink(output_dir_id_no_msi_pole, recursive = TRUE)
syn_exposures_id_no_msi_pole <-
  syn_tumors_id_no_msi_pole$ground.truth.exposures

# Generate MSI-H synthetic tumors
synthetic_tumors_id_msi <-
  generate_subtype_syn_tumors(
    seed = seed,
    dir = output_dir_id_msi,
    cancer_types = cancer_types_msi,
    samples_per_caner_type = num_samples_msi_scaled,
    input_sigs = input_sigs_id,
    real_exposure = real_exposures_id_msi,
    distribution = distribution,
    sample_prefix_name = sample_prefix_name,
    tumor_marker_name = "MSI-H",
    sig_params = sig_params_id_nine_types
  )
unlink(output_dir_id_msi, recursive = TRUE)
syn_exposures_id_msi <-
  synthetic_tumors_id_msi$ground.truth.exposures

# Generate POLE-mutated synthetic tumors
synthetic_tumors_id_pole <-
  generate_subtype_syn_tumors(
    seed = seed,
    dir = output_dir_id_pole,
    cancer_types = cancer_types_pole,
    samples_per_caner_type = num_samples_pole_scaled,
    input_sigs = input_sigs_id,
    real_exposure = real_exposures_id_pole,
    distribution = distribution,
    sample_prefix_name = sample_prefix_name,
    tumor_marker_name = "POLE",
    sig_params = sig_params_id_nine_types
  )
unlink(output_dir_id_pole, recursive = TRUE)
syn_exposures_id_pole <-
  synthetic_tumors_id_pole$ground.truth.exposures

# Combine the non MSI-H non POLE, MSI-H and POLE synthetic exposures in each
# cancer type
synthetic_exposures_id <-
  combine_exposure(
    syn_exposures_id_no_msi_pole,
    syn_exposures_id_pole,
    syn_exposures_id_msi
  )

# Generate the combined synthetic tumors
write_sig_params(
  dir = output_dir_id,
  real_exposure = real_exposures_id,
  synthetic_exposure = synthetic_exposures_id,
  cancer_types = cancer_types,
  distribution = distribution,
  sig_params = sig_params_id_nine_types,
  sample_prefix_name = sample_prefix_name,
  mutation_type = mutation_type
)

catalog.info <- SynSigGen:::CreateSynCatalogs(
  signatures = input_sigs_id,
  exposures  = synthetic_exposures_id)
ICAMS::WriteCatalog(
  catalog.info$ground.truth.catalog, 
  file = file.path(output_dir_id, "ground.truth.syn.catalog.csv"))
mSigTools::write_exposure(
  exposure = synthetic_exposures_id, 
  file     = file.path(output_dir_id, "ground.truth.syn.exposures.csv"))

# Add noise to the synthetic tumors
id_noisy_tumors_size_selected <-
  SynSigGen::GenerateNoisyTumors(
    seed = seed,
    dir = output_dir_id_nb_size_selected,
    input.exposure = synthetic_exposures_id,
    signatures = input_sigs_id,
    n.binom.size = 4
  )

noisy_exposures_size_selected_id <- 
  id_noisy_tumors_size_selected$exposures

get_sig_universes(
  exposures = real_exposures_id, 
  filename = file.path(output_dir_id_nb_size_selected,
                       "ground.truth.sig.universe.csv"),
  sigs = input_sigs_id,
  sig_file = file.path(output_dir_id_nb_size_selected,
                       "ground.truth.sigs.csv"))

#################################################################
##                   Plot data distributions                   ##
#################################################################

data_distribution_file <-
  "./synthetic_data/ID/intermed_results/ID_syn_data_distributions.pdf"
grDevices::pdf(
  file = data_distribution_file,
  width = 11.6929, height = 8.2677, onefile = TRUE
)
par(mfrow = c(3, 3))
plot_exposure_distribution(
  real_exposure = real_exposures_id,
  synthetic_exposure = synthetic_exposures_id,
  noisy_exposure = noisy_exposures_size_selected_id,
  size = 4,
  distribution = distribution,
  sig_params = sig_params_id_nine_types,
  sample_prefix_name = sample_prefix_name
)
grDevices::dev.off()

#################################################################
##               Put data in standard places / file names      ##
#################################################################

source("synthetic_data/data_gen_code/data_gen_rename.R")

old_dataset_names <- basename(output_dir_id_nb_size_selected)

data_gen_rename(dataset = "ID",
                old_dataset_name = old_dataset_names, 
                top_folder_name = "synthetic_data")
