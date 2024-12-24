# What are the cosine similarities of the reconstructed PCAWG exposures
# for stomach ?

if(!grepl("attribution_paper_code", basename(getwd()))) {
  stop("Run this script from the top level directory")
}
rm(list = ls())

source("analysis/code/reoptimize_similarities_from_exposures.R")
library(dplyr)
library(tibble)
library(ggplot2)
library(huxtable)
library(tidyverse)

if (FALSE) { # Not used
  synthetic_stomach_reoptimized = 
    reoptimize_similarities_from_exposures(
      synthetic_stomach_spectra(),
      synthetic_stomach_exposures())
  save(
    synthetic_stomach_reoptimized, 
    file = "output_for_paper/synthetic_stomach_reoptimized.rdata")
}

if (FALSE) { # Run the code in this next block to generate syn_all_reoptimized.rdata
syn_all_reoptimized = reoptimize_similarities_from_exposures(
  get_all_input("SBS")$all_spectra,
  get_ground_truth_exposure("SBS")
)
save(syn_all_reoptimized,
     file =
       "output_for_paper/syn_all_reoptimized.rdata")
}
#######

load("output_for_paper/syn_all_reoptimized.rdata")
dplyr::select(syn_all_reoptimized, 
              sample_id, 
              cosine_similarity = cossim_for_NNLS) %>%
  rename(`Sample ID` = sample_id, `Cosine similarity` = cosine_similarity) %>%
  huxtable::as_huxtable() %>%
  set_wrap(TRUE) %>%
  set_align("center") %>%
  set_valign("middle") %>%
  set_escape_contents(1:nrow(.), 1, TRUE) %>%
  set_number_format(2:nrow(.), 2, "%5.4f") %>%
  insert_row(
    paste("Table S2: For all synthetic data, the cosine similarity",
           "of the reconstruction using the ground-truth attribution",
           "signatures to the synthetic spectrum"),
    fill = "",
    after = 0) %>%
  set_row_height(1, 3 / nrow(.)) %>%
  set_row_height(2, 1.5 / nrow(.)) %>%
  set_wrap(row = 1, value = TRUE) %>%
  merge_cells(1, 1:2) %>%
  set_header_rows(1, TRUE) %>%
  style_header_rows(bold = TRUE) %>%
  quick_xlsx(
    file = "output_for_paper/table_s2_cossim_for_synthetic.xlsx")
