if(!grepl("attribution_paper_code", basename(getwd()))) {
  stop("Run this script from the top level directory")
}
rm(list = ls())
message("ls() = ", ls())
source("analysis/code/powerset_analysis.R")
library(huxtable)
library(openxlsx)


all_stomach_spectra <- pcawg_stomach_spectra()

# bigfile = "~/bigdata/powerset_real_spectra.rdata"
bigfile = "output_for_paper/powerset_real_spectra.rdata"
message("Loading existing file ", bigfile)
load(bigfile)
stopifnot(length(powerset_real_spectra) == 75)

# 
if (is.null(names(powerset_real_spectra))) {
  names(powerset_real_spectra) = colnames(all_stomach_spectra)
}
#

all_exp <- lapply(powerset_real_spectra, `[[`, "exposures")
num_recon <- unlist(lapply(all_exp, number_reconstructions))

sum(num_recon > 0)
# 65
median(num_recon)
# 11189

median(num_recon[num_recon > 0])
# 257

figure_2_histogram <- function(file, num_of_recon) {
  cairo_pdf(file = file)
  par(mfrow = c(2, 1), mai = c(1.02, 1.25, 0.82, 0.42))
  hist(num_of_recon,
       breaks = seq(from = 0, to = 68000, by = 1000),
       xlab = paste0("\nNumber of attributions based on distinct sets of signatures\n",
                     "that yield cosine similarity to spectrum > 0.969"),
       ylab = "Number of stomach\nadenocarcinoma spectra",
       cex.lab = 0.9,
       cex.axis = 0.9,
       main = ""
  )
  dev.off()
}

figure_2_histogram(
  file = "output_for_paper/fig_2A.pdf",
  num_of_recon = num_recon[num_recon > 0]
)

check_one_sample_for_no_low_fraction = function(one_sample_name, fraction) {
  one_sample = powerset_real_spectra[[one_sample_name]]
  result1 = unlist(
    lapply(X = one_sample$exposures, 
           filter_exposure_no_low_fraction_contribution,
           fraction = fraction))
  # browser()
  return(list(sample_name = one_sample_name, num_recon = sum(result1)))
}

counts_at_fraction_threshold_0.06 = 
  mSigAct:::ListOfList2Tibble(
    lapply(names(powerset_real_spectra),
             check_one_sample_for_no_low_fraction, fraction = 0.06))
sum(counts_at_fraction_threshold_0.06$num_recon > 1)
# [1] 47
mean(counts_at_fraction_threshold_0.06$num_recon)
# [1] 38.32

counts_at_fraction_threshold_0.03 = 
  mSigAct:::ListOfList2Tibble(
    lapply(names(powerset_real_spectra),
           check_one_sample_for_no_low_fraction, fraction = 0.03))
sum(counts_at_fraction_threshold_0.03$num_recon > 1)
# 64
mean(counts_at_fraction_threshold_0.03$num_recon)
# 447.4
sum(counts_at_fraction_threshold_0.03$num_recon == 0)
# 11

# Generate main text figure showing manually selected example reconstructions
# and generate corresponding draft supplementary table
example_sample_name <- "Stomach-AdenoCA::SP85251"

dplyr::filter(counts_at_fraction_threshold_0.03, 
              sample_name == example_sample_name)
# 120

dplyr::filter(counts_at_fraction_threshold_0.06, 
              sample_name == example_sample_name)
# 8


example_exposures = powerset_real_spectra[[example_sample_name]]
 # check_one_sample_for_no_low_fraction(example_sample_name, 0.06)
attribution_indices_to_show =
  which(unlist(lapply(example_exposures$exposures, filter_exposure_no_low_fraction_contribution, 0.03)))

attributions_to_show = example_exposures$exposures[attribution_indices_to_show]
reconstructions_to_show = example_exposures$reconstructions[ , attribution_indices_to_show]

lapply(attributions_to_show, `[[`, "q_exp_matrix") %>%
  lapply(data.frame) %>%
  lapply(function(ex) { 
    colnames(ex) = paste(rownames(ex), collapse = "_")
    return(cbind(signature_id = rownames(ex), ex))}
  ) -> column_list
                                                   
purrr::reduce(.x = column_list, .f = dplyr::full_join) %>%
  dplyr::slice(gtools::mixedorder(signature_id)) ->
  table_v0
table_v1 = table_v0[ , -1]

unlist(lapply(attributions_to_show, `[[`, 4)) -> cosines
names(cosines) = colnames(table_v1)

rbind(cosines, table_v1) %>%
  `rownames<-`(c(
  paste("cosine similarity to spectrum of", example_sample_name),
  table_v0[ , 1])) %>%
  huxtable::as_huxtable(add_rownames = TRUE, add_colnames = FALSE) %>%
  insert_row(
    paste("Table S4: For Stomach-AdenoCA from Alexandrov et al, 2020,",
          "the cosine similarity",
          "of the reconstruction using the attributions in that",
          "paper to the actual spectrum"),
    fill = "",
    after = 0) %>%
  set_wrap(row = 1, value = TRUE) %>%
  merge_cells(1, 1:ncol(.)) %>%
  set_escape_contents(1, 1:ncol(.), TRUE) %>%
  set_header_rows(1, TRUE) %>%
  style_header_rows(bold = TRUE) %>%
  set_contents(2, 1, "cosine similarity between reconstruction and spectrum") %>%
  set_wrap(2, 1, TRUE) %>%
  set_align(2:nrow(.), 1:ncol(.), "center") %>%
  set_valign("middle") %>%
  set_number_format(2, 2:ncol(.), "%5.4f") %>%
  huxtable::quick_xlsx(
    file = "output_for_paper/table_s4_pcawg_stomach_powerset.xlsx")

plotrange = 1:5  
plot_reconstruction_pdf(
  file = "output_for_paper/fig_2B-G.pdf",
  spectrum = all_stomach_spectra[ , example_sample_name, drop = FALSE],
  reconstructions = reconstructions_to_show[ , plotrange],
  exposures = attributions_to_show[plotrange],
  xlabels = FALSE

)

