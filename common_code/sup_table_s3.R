if(basename(getwd()) != "sig_attribution_paper_code") {
  stop("Run this script from the top level directory")
}

rm(list = ls())

source("analysis/code/analysis_utils.R")

all_stomach <- mSigAct::ExposureProportions(
  mutation.type = "SBS96",
  cancer.type = "Stomach-AdenoCA"
)
selected_sig_names <- names(all_stomach[all_stomach > min(all_stomach)])

all_syn <-
  ICAMS::ReadCatalog("synthetic_data/SBS/ground.truth.syn.catalog.csv")
all_syn_exp <-
  data.table::fread("synthetic_data/SBS/ground.truth.syn.exposures.csv")

# These are the signatures removed before
# exploring the power set: "SBS9"  "SBS20" "SBS28" "SBS43" "SBS51" "SBS58"
# Among these, the synthetic data contained only SBS20, SBS28
keep <- c(1, grep("Stomach-AdenoCA", (colnames(all_syn_exp))))
stom.exp <- all_syn_exp[, ..keep]
st <- t(stom.exp[, 2:101])
colnames(st) <- dplyr::pull(stom.exp, V1)

# None of the synthetic stomach spectra had SBS28
stopifnot(all(st[, "SBS28"] == 0))

# Remove synthetic stomach spectra with SBS20
st2 <- st[st[, "SBS20"] == 0, ]
stopifnot(dim(st2) == c(95, 36))
# st2 now contains the synthetic Stomach-AdenoCA that do not have SBS20
# The rows of st2 are samples and the columns are signatures

if (!exists("sup_table_s3_results")) {
  if (file.exists("common_code/sup_table_s3_results.Rdata")) {
    load("common_code/sup_table_s3_results.Rdata")
  } else {
    sup_table_s3_results <-
      get_alternative_attribution(
        sample_names = rownames(st2),
        spectra = all_syn,
        selected_sig_names = selected_sig_names
      )
    save(sup_table_s3_results, file = "common_code/sup_table_s3_results.Rdata")
  }
}
table_results <- get_alter_attribution_table(sup_table_s3_results)

nrow(table_results) # 95
sum(table_results$num_all_better > 0) # 93
median(table_results$num_all_better) # 15

output_home <- "output_for_paper/"
writexl::write_xlsx(
  x = table_results,
  path = file.path(output_home, "sup_table_s3.xlsx")
)
