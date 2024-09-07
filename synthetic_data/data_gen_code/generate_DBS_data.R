# Run this script from the top level directory

source("synthetic_data/data_gen_code/data_gen_utils.R")

library(cosmicsig)
library(PCAWG7)
library(SynSigGen)

##################################################################
##                      Data preprocessing                      ##
##################################################################
# Get the real exposures from PCAWG assignments
real_exposures_dbs78 <- PCAWG7::exposure$PCAWG$DBS78
pcawg_dbs78_catalog <- PCAWG7::spectra$PCAWG$DBS78

# Only select samples that belong to the selected cancer types
cancer_types <- c(
  "Breast-AdenoCA", "ColoRect-AdenoCA", "Eso-AdenoCA",
  "Kidney-RCC", "Liver-HCC", "Lung-AdenoCA",
  "Ovary-AdenoCA", "Skin-Melanoma", "Stomach-AdenoCA"
)
indices_nine_types <- unlist(sapply(cancer_types, FUN = function(x) {
  grep(x, colnames(real_exposures_dbs78))
}))
real_exposures_dbs78 <- real_exposures_dbs78[, indices_nine_types]

# Exclude samples which have mutations less than 100
samples_less_than_100 <- names(which(colSums(pcawg_dbs78_catalog) < 100))
indices_less_than_100 <-
  which(colnames(real_exposures_dbs78) %in% samples_less_than_100)
real_exposures_dbs78 <-
  real_exposures_dbs78[, -indices_less_than_100, drop = FALSE]

real_exposures_dbs78 <- remove_zero_activity_sigs(real_exposures_dbs78)

##################################################################
##   Calculate number of synthetic tumors in each cancer type   ##
##################################################################

msi_sample_indices_nine_types <-
  unlist(sapply(pcawg_msi_tumor_ids, FUN = function(x) {
    grep(x, colnames(real_exposures_dbs78))
  }))
msi_sample_ids <- names(msi_sample_indices_nine_types)
length(msi_sample_ids) # 15

pole_sample_indices_nine_types <-
  unlist(sapply(pcawg_pole_tumor_ids, FUN = function(x) {
    grep(x, colnames(real_exposures_dbs78))
  }))
pole_sample_ids <- names(pole_sample_indices_nine_types)
length(pole_sample_ids) # 6

real_exposures_dbs78_no_msi_pole <-
  real_exposures_dbs78[, -c(
    msi_sample_indices_nine_types,
    pole_sample_indices_nine_types
  ), drop = FALSE]
real_exposures_dbs78_no_msi_pole <-
  remove_zero_activity_sigs(real_exposures_dbs78_no_msi_pole)
real_exposures_dbs78_msi <-
  real_exposures_dbs78[, msi_sample_indices_nine_types, drop = FALSE]
real_exposures_dbs78_msi <- remove_zero_activity_sigs(real_exposures_dbs78_msi)
real_exposures_dbs78_pole <-
  real_exposures_dbs78[, pole_sample_indices_nine_types, drop = FALSE]
real_exposures_dbs78_pole <- remove_zero_activity_sigs(real_exposures_dbs78_pole)

num_samples_total <- calculate_num_samples(real_exposures_dbs78)
num_samples_msi <- calculate_num_samples(real_exposures_dbs78_msi)
cancer_types_msi <- names(num_samples_msi)
num_samples_pole <- calculate_num_samples(real_exposures_dbs78_pole)
cancer_types_pole <- names(num_samples_pole)
num_samples_no_msi_pole <- calculate_num_samples(real_exposures_dbs78_no_msi_pole)

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

# Make sure the total number of synthetic tumors is 900
sum(num_samples_msi_scaled) + sum(num_samples_pole_scaled) +
  sum(num_samples_no_msi_pole_scaled)

##################################################################
##                 Generation of synthetic data                 ##
##################################################################

output_dir_dbs78_no_msi_pole <- "./synthetic_data/DBS/intermed_results/DBS78.syn.exposures.no.msi"
output_dir_dbs78_msi <- "./synthetic_data/DBS/intermed_results/DBS78.syn.exposures.msi"
output_dir_dbs78_pole <- "./synthetic_data/DBS/intermed_results/DBS78.syn.exposures.pole"
output_dir_dbs78 <- "./synthetic_data/DBS/intermed_results/DBS78.syn.exposures.no.noise"
output_dir_dbs78_nb_size_selected <-
  "./synthetic_data/DBS/intermed_results/DBS78.syn.exposures.noisy.neg.binom.size.selected"

distribution <- "neg.binom"
sample_prefix_name <- ""
mutation_type <- "DBS78"
seed <- 658220
input_sigs_dbs78 <- cosmicsig::COSMIC_v3.4$signature$GRCh37$DBS78

sig_params_dbs78_nine_types <-
  SynSigGen:::GetSynSigParamsFromExposures(
    exposures = real_exposures_dbs78,
    distribution = distribution,
    sig.params = SynSigGen::signature.params$DBS78
  )

# Generate synthetic tumors that are not MSI-H or POLE-mutated
synthetic_tumors_dbs78_no_msi_pole <-
  SynSigGen::GenerateSyntheticTumors(
    seed = seed,
    dir = output_dir_dbs78_no_msi_pole,
    cancer.types = cancer_types,
    samples.per.cancer.type = num_samples_no_msi_pole_scaled,
    input.sigs = input_sigs_dbs78,
    real.exposures = real_exposures_dbs78_no_msi_pole,
    distribution = distribution,
    sample.prefix.name = sample_prefix_name,
    sig.params = sig_params_dbs78_nine_types
  )
unlink(output_dir_dbs78_no_msi_pole, recursive = TRUE)
syn_exposures_dbs78_no_msi_pole <-
  synthetic_tumors_dbs78_no_msi_pole$ground.truth.exposures

# Generate MSI-H synthetic tumors
synthetic_tumors_dbs78_msi <-
  generate_subtype_syn_tumors(
    seed = seed,
    dir = output_dir_dbs78_msi,
    cancer_types = cancer_types_msi,
    samples_per_caner_type = num_samples_msi_scaled,
    input_sigs = input_sigs_dbs78,
    real_exposure = real_exposures_dbs78_msi,
    distribution = distribution,
    sample_prefix_name = sample_prefix_name,
    tumor_marker_name = "MSI-H",
    sig_params = sig_params_dbs78_nine_types
  )
unlink(output_dir_dbs78_msi, recursive = TRUE)
syn_exposures_dbs78_msi <-
  synthetic_tumors_dbs78_msi$ground.truth.exposures

# Generate POLE-mutated synthetic tumors
synthetic_tumors_dbs78_pole <-
  generate_subtype_syn_tumors(
    seed = seed,
    dir = output_dir_dbs78_pole,
    cancer_types = cancer_types_pole,
    samples_per_caner_type = num_samples_pole_scaled,
    input_sigs = input_sigs_dbs78,
    real_exposure = real_exposures_dbs78_pole,
    distribution = distribution,
    sample_prefix_name = sample_prefix_name,
    tumor_marker_name = "POLE",
    sig_params = sig_params_dbs78_nine_types
  )
unlink(output_dir_dbs78_pole, recursive = TRUE)
syn_exposures_dbs78_pole <-
  synthetic_tumors_dbs78_pole$ground.truth.exposures

# Combine the non MSI-H non POLE, MSI-H and POLE synthetic exposures in each
# cancer type
synthetic_exposures_dbs78 <-
  combine_exposure(
    syn_exposures_dbs78_no_msi_pole,
    syn_exposures_dbs78_pole,
    syn_exposures_dbs78_msi
  )

# Generate the combined synthetic tumors
write_sig_params(
  dir = output_dir_dbs78,
  real_exposure = real_exposures_dbs78,
  synthetic_exposure = synthetic_exposures_dbs78,
  cancer_types = cancer_types,
  distribution = distribution,
  sig_params = sig_params_dbs78_nine_types,
  sample_prefix_name = sample_prefix_name,
  mutation_type = mutation_type
)

catalog.info <- SynSigGen:::CreateSynCatalogs(
  signatures = input_sigs_dbs78,
  exposures  = synthetic_exposures_dbs78)
ICAMS::WriteCatalog(
  catalog.info$ground.truth.catalog, 
  file = file.path(output_dir_dbs78, "ground.truth.syn.catalog.csv"))
mSigTools::write_exposure(
  exposure = synthetic_exposures_dbs78, 
  file     = file.path(output_dir_dbs78, "ground.truth.syn.exposures.csv"))

# Add noise to the synthetic tumors
dbs78_noisy_tumors_size_selected <-
  SynSigGen::GenerateNoisyTumors(
    seed = seed,
    dir = output_dir_dbs78_nb_size_selected,
    input.exposure = synthetic_exposures_dbs78,
    signatures = input_sigs_dbs78,
    n.binom.size = 1
  )

noisy_exposures_size_selected_dbs78 <- 
  dbs78_noisy_tumors_size_selected$exposures

# duplicate call because downstream copying code expects 
# files in both locations
get_sig_universes(
  exposures = real_exposures_dbs78, 
  filename = file.path(output_dir_dbs78_nb_size_selected,
                       "ground.truth.sig.universe.csv"),
  sigs = input_sigs_dbs78,
  sig_file = file.path(output_dir_dbs78_nb_size_selected,
                       "ground.truth.sigs.csv"))

#################################################################
##                   Plot data distributions                   ##
#################################################################

data_distribution_file <-
  "./synthetic_data/DBS/intermed_results/DBS_syn_data_distribution.pdf"
grDevices::pdf(
  file = data_distribution_file,
  width = 11.6929, height = 8.2677, onefile = TRUE
)
par(mfrow = c(3, 3))
plot_exposure_distribution(
  real_exposure = real_exposures_dbs78,
  synthetic_exposure = synthetic_exposures_dbs78,
  noisy_exposure = noisy_exposures_size_selected_dbs78,
  size = 1,
  distribution = distribution,
  sig_params = sig_params_dbs78_nine_types,
  sample_prefix_name = sample_prefix_name
)
grDevices::dev.off()

#################################################################
##               Put data in standard places / file names      ##
#################################################################

source("synthetic_data/data_gen_code/data_gen_rename.R")

old_dataset_names <- basename(output_dir_dbs78_nb_size_selected)

data_gen_rename(dataset = "DBS",
                old_dataset_name = old_dataset_names, 
                top_folder_name = "synthetic_data")
