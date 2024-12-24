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


if (FALSE) {
  pcawg_stomach_reoptimized = 
    reoptimize_similarities_from_exposures(
      pcawg_stomach_spectra(),
      pcawg_stomach_exposures())
  save(
    pcawg_stomach_reoptimized, 
    file = "output_for_paper/pcawg_stomach_reoptimized.rdata")
}

pcawg_all_reoptimized = 
  reoptimize_similarities_from_exposures(
    PCAWG7::spectra$PCAWG$SBS96,
    PCAWG7::exposure$PCAWG$SBS96)
save(
  pcawg_all_reoptimized, 
  file = "output_for_paper/pcawg_all_reoptimized.rdata")

median(pcawg_all_reoptimized$cossim_for_NNLS)
# [1] 0.9688108

######

load("output_for_paper/pcawg_all_reoptimized.rdata")
dplyr::select(pcawg_all_reoptimized, 
              sample_id, 
              cosine_similarity = cossim_for_NNLS) %>%
  rename(`Sample ID` = sample_id, `Cosine similarity` = cosine_similarity) %>%
  huxtable::as_huxtable() %>%
  set_escape_contents(1:nrow(.), 1, TRUE) %>%
  set_number_format(2:nrow(.), 2, "%5.4f") %>%
  insert_row(
    paste("Table S1: For all data from Alexandrov et al, 2020,",
          "the cosine similarity",
          "of the reconstruction using the attributions",
          "in that paper to the actual spectrum"),
    fill = "",
    after = 0) %>%
  set_wrap(row = 1, value = TRUE) %>%
  merge_across(1, 1:2) %>%
  set_header_rows(1, TRUE) %>%
  style_header_rows(bold = TRUE) %>%
  set_row_height(1, 3 / nrow(.)) %>%
  set_row_height(2, 1.5 / nrow(.)) %>%
  set_align(2:nrow(.), 1:ncol(.), "center") %>%
  set_valign("middle") %>%
  quick_xlsx(
    file = "output_for_paper/table_s1_cossim_for_pcawg.xlsx")



######


if (FALSE) {
  pcawg_stomach_reoptimized %>% filter(logLH > -10000) ->
    p_filtered
  
  
  ggplot(p_filtered, aes(x = logLH)) +
    geom_histogram(binwidth = 10, color = "black", fill = "blue") +
    scale_x_continuous(breaks = seq(-10000, 0, by = 500)) +
    labs(title = "Histogram of logLH (Filtered > -10000)",
         x = "logLH",
         y = "Count") +
    theme_minimal()
  
  library(rlang)  # Make sure to load rlang for sym()
  
  similarity_scatter <- function(tb, xcol, ycol, title_info) {
    # Convert column names (strings) to symbols
    xcol_sym <- sym(xcol)
    ycol_sym <- sym(ycol)
    
    ggplot(tb, aes(x = !!xcol_sym, y = !!ycol_sym)) +
      geom_point(color = "blue", size = 2) +
      labs(
        title = paste("Scatter Plot of cosine similarity vs logLH for", title_info, "Stomach-AdenoCA"),
        x = xcol,  # Use the string itself for labeling
        y = ycol   # Use the string itself for labeling
      ) +
      theme_minimal()
  }
  
  # Example call
  similarity_scatter(p_filtered, "logLH", "cossim_for_NNLS", "PCAWG")
}


