library(ggplot2)
library(ggbeeswarm)
library(gridExtra)
library(dplyr)
library(data.table)
library(huxtable)

ranks_by_tool = function(mutation_type, mean_best_worst) {
  ranktable = fread(
    paste0("output_for_paper/tools_ranked_in_each_cancer_type_",
    mutation_type,
    ".csv"))

  ranktable = group_by(ranktable, Tool)
  
  if (mean_best_worst == "mean") {
    myfn = mean
    newcolname = "Mean rank"
    partfilename = "_mean"
  } else if (mean_best_worst == "best") {
    myfn = min
    newcolname = "Best rank"
    partfilename = "_best"
  } else if (mean_best_worst == "worst") {
    myfn = max
    newcolname = "Worst rank"
    partfilename = "_worst"
  } else if (mean_best_worst == "sd") {
    myfn = sd
    newcolname = "Standard deviation of rank"
    partfilename = "_sd"
  } else {
    stop("unrecognized measure: ", mean_best_worst)
  }

  means <- ranktable %>%
    group_by(Tool) %>%
    summarize(target = myfn(`Rank within cancer type`)) %>%
    ungroup() %>%
    arrange(target)
  
  
  means %>%
    mutate(`Cancer type` = "", rank = NULL) %>%
    relocate(`Cancer type`, Tool, target) ->
    table2

    colnames(table2) = c("", "Tool", newcolname)

    table2 %>%
      huxtable::as_huxtable() %>%
      huxtable::set_align("centre") %>%
      huxtable::set_number_format(NA) %>%
      huxtable::insert_row(
        paste("Summary,", mean_best_worst, "ranks for",
              mutation_type, "signatures"),
        fill = "") %>%
      huxtable::merge_across(1, 1:3) %>%
      huxtable::set_header_rows(1:2, TRUE) %>%
      huxtable::style_header_rows(bold = TRUE) %>%
      quick_xlsx(
        file = 
          file.path("output_for_paper",
                    paste0("add_this_manually_to_sup_table_for_",
                           mutation_type,
                           partfilename,
                           "_rank.xlsx")))
  
}

for (mutation_type in c("SBS", "DBS", "ID")) {
  for (rank_type in c("mean", "best", "worst", "sd")) {
    ranks_by_tool(mutation_type, rank_type)
  }
}

