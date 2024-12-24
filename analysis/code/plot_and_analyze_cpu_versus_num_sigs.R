# Source this file to analyze and plot the relationship between
# CPU seconds and the number of signatures considered for attribution
# in a given cancer type.

rm(list = ls())

source("analysis/code/common_utils.R")
source("analysis/code/get_all_input.R")
library(ggplot2)
library(dplyr)
library(purrr)
library(broom)
library(ggforce)
library(huxtable)

set.seed(1066)
cpu_versus_num_sigs_table = 23

cpu_versus_num_sigs = function(mutation_type) {
  all_inputs = get_all_input(mutation_type)

  num_sigs = lapply(all_inputs$signature_universes, length)
  num_sigs_by_cancer_type = tibble(cancer_type = names(num_sigs), num_sigs = unlist(num_sigs))
  # browser()
  
  base_table_num = 20
  prefix = "output_for_paper/table_s"
  suffix = "_cpu_by_cancer_type_"
  if (mutation_type == "SBS") {
    timings = data.table::fread(
      paste0(prefix, base_table_num, suffix,"SBS.csv"))
  } else if (mutation_type == "DBS") {
    timings = data.table::fread(
      paste0(prefix, base_table_num + 1, suffix, "DBS.csv"))
  } else if (mutation_type == "ID") {
    timings = data.table::fread(
      paste0(prefix, base_table_num + 2, suffix, "ID.csv"))
  }
  df = dplyr::full_join(num_sigs_by_cancer_type, timings)
  mutate(df, tool_name = pretty_tool_names(tool_name)) %>%
    filter(cancer_type != "Total all cancer types") -> df
  
  # browser()
  if (FALSE) {
  for (tool in unique(df$tool_name)) {
    tmpdf = filter(df, tool_name == tool)
    plot(cpu_seconds ~ num_sigs, data = tmpdf, main = tool)
    
  }
  }
  
  correlation_results <- df %>%
    group_by(tool_name) %>%
    tidyr::nest() %>%
    mutate(
      correlation_test = map(data, ~ cor.test(jitter(.x$num_sigs), jitter(.x$cpu_seconds), method = "spearman")),
      tidy_results = map(correlation_test, broom::tidy)
    ) %>%
    tidyr::unnest(tidy_results) %>%
    select(tool_name, estimate, p.value)

  correlation_results$FDR = p.adjust(correlation_results$p.value, method = "BH")
  correlation_results$Mutation_type = mutation_type
  print(correlation_results)
  
  nrow <- 5
  ncol <- 3
  n_plots_per_page <- nrow * ncol
  n_pages <- ceiling(length(unique(df$tool_name)) / n_plots_per_page)
  
  sup_fig_nums = c(SBS = 8, DBS = 9, ID = 10)
  
  pdf(file = paste0("output_for_paper/sup_fig_S",
                    sup_fig_nums[mutation_type], "_", 
                    mutation_type, 
                    "_cpu_time_vs_sig_count.pdf"))
  # Loop over pages and print plots for each page
  for (i in 1:n_pages) {
    plot_object <- ggplot(df, aes(x = num_sigs, y = cpu_seconds)) +
      geom_jitter() +
      facet_wrap_paginate(~ tool_name, nrow = nrow, ncol = ncol, page = i, scales = "free", shrink = FALSE) +  # Pagination
      labs(x = paste("Number of", mutation_type, "signatures in cancer type"),
           y = paste("CPU Time (seconds) to attribute", mutation_type, "signatures")) +
      theme_minimal() +
      scale_x_continuous(breaks = seq(min(df$num_sigs, na.rm = TRUE),
                                      max(df$num_sigs, na.rm = TRUE),
                                      by = 1)) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate labels 90 degrees
      )
    
    print(plot_object)  # Print each page of plots
  }
  dev.off()
  return(correlation_results)
}

correlations = lapply(c("SBS", "DBS", "ID"), cpu_versus_num_sigs)
do.call(rbind, correlations) -> correlations
lowfdr = which(correlations$FDR < 0.1)

correlations %>%
  ungroup() %>%
  mutate(
    `Mutation type` = Mutation_type,
    Estimate = estimate, P = p.value, Tool = tool_name, Estimate = estimate) %>%
  select(`Mutation type`, Tool, P, FDR) %>%
  huxtable::as_huxtable() %>%
  set_align("center") %>%
  set_wrap(TRUE) %>%
  insert_row(
    paste("P value by Spearman rank correlation,",
          "FDR by Benjamini-Hochberg method within each mutation type"),
    fill = ""
  ) %>%
  insert_row(
    paste0("Table S", cpu_versus_num_sigs_table,
           ": CPU seconds as a function of ",
          "the number of mutational signatures in a cancer type."),
          fill = ""
  ) %>%
  set_bold(lowfdr + 3, 1:ncol(.)) %>%
  set_header_rows(c(1, 3), TRUE) %>%
  style_header_rows(bold = TRUE) %>%
  merge_across(1:2, 1:ncol(.)) %>%
  set_row_height(row = 1:2, value = 1.8 / nrow(.)) %>%
  set_row_height(row = 3, value = 1.4 / nrow(.)) %>%
  set_valign("middle") %>%
  quick_xlsx(
    file = paste0("output_for_paper/table_s",
                  cpu_versus_num_sigs_table, 
                  "_cpu_versus_number_of_sigs.xlsx"))
