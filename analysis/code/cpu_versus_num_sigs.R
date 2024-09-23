source("analysis/code/get_all_input.R")
library(ggplot2)
library(dplyr)
library(purrr)
library(broom)
library(ggforce)

cpu_versus_num_sigs = function(mutation_type) {
  all_inputs = get_all_input(mutation_type)

  num_sigs = lapply(all_inputs$signature_universes, length)
  num_sigs_by_cancer_type = tibble(cancer_type = names(num_sigs), num_sigs = unlist(num_sigs))
  browser()
  timings = data.table::fread("analysis/summary/SBS/cpu_seconds_by_cancer_type_SBS.csv")
  df = dplyr::full_join(num_sigs_by_cancer_type, timings)
  df = mutate(df, base_tool = gsub("_0.*", "", tool_name)) 
  
  correlation_results <- df %>%
    group_by(base_tool) %>%
    nest() %>%
    mutate(
      correlation_test = map(data, ~ cor.test(.x$num_sigs, .x$cpu_seconds)),
      tidy_results = map(correlation_test, tidy)
    ) %>%
    unnest(tidy_results) %>%
    select(base_tool, estimate, p.value)

  correlation_results$adj.p = p.adjust(correlation_results$p.value, method = "BH")
  print(correlation_results)
  
  nrow <- 3
  ncol <- 2
  n_plots_per_page <- nrow * ncol
  n_pages <- ceiling(length(unique(df$base_tool)) / n_plots_per_page)
  
  # These labels are only for SBS
  custom_x_labels <- c(
    "Breast-AdenoCA" = 13,
    "Eso-AdenoCA" = 10,
    "Kidney-RCC" = 14,
    "Stomach-AdenoCA" = 16,
    "Liver-HCC" = 20
  )
  
  # Loop over pages and print plots for each page
  for (i in 1:n_pages) {
    plot_object <- ggplot(df, aes(x = num_sigs, y = cpu_seconds)) +
      geom_jitter() +
      facet_wrap_paginate(~ base_tool, nrow = nrow, ncol = ncol, page = i, scales = "free_y") +  # Pagination
      labs(x = "Number of Signatures", y = "CPU Time (seconds)") +
      theme_minimal() +
      scale_x_continuous(
        breaks = custom_x_labels,  # Breaks at specific num_sigs values
        labels = names(custom_x_labels)  # Corresponding labels
      ) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)  # Rotate labels 90 degrees
      )
    
    print(plot_object)  # Print each page of plots
  }

}

cpu_versus_num_sigs("DBS")