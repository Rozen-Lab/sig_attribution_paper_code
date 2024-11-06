source("analysis/code/get_all_input.R")
library(ggplot2)
library(dplyr)
library(purrr)
library(broom)
library(ggforce)

set.seed(1066)

score_versus_num_sigs = function(mutation_type) {
  all_inputs = get_all_input(mutation_type)

  num_sigs = lapply(all_inputs$signature_universes, length)
  num_sigs_by_cancer_type = tibble(cancer_type = names(num_sigs), num_sigs = unlist(num_sigs))
  # browser()
  
  
  scores = data.table::fread(paste0("output_for_paper/stats_each_sample_", mutation_type, ".csv"))
  
  scores2 = dplyr::mutate(scores, cancer_type = gsub("::.*", "", Sample.ID))
  

  dplyr::full_join(num_sigs_by_cancer_type, scores2) %>%
    filter(Tool %in% global_raw_tools_to_plot) %>%
    mutate(tool_name = pretty_tool_names(Tool)) %>%
    filter(cancer_type != "Total all cancer types") -> df
  
  
  pdf(paste0("output_for_paper/tmp_score_vs_signum_", mutation_type, ".pdf"))
  par(mfrow = c(2, 1))
  tool_p = 
    lapply(unique(df$tool_name),
           function(tooln) {
             tmpdf = filter(df, tool_name == tooln)
             message(tooln)
             # browser()
             mm = lm(Combined ~ num_sigs, data = tmpdf)
             boxplot(Combined ~ num_sigs, data = tmpdf, main = tooln)
             
             list(`Mutation type` = mutation_type,
                  Tool = tooln, 
                  `Estimated coefficient` = summary(mm)$coefficients["num_sigs", 1],
                  P = summary(mm)$coefficients["num_sigs", 4])
           })
  dev.off()
  
  res = mSigAct:::ListOfList2Tibble(tool_p)
  res$FDR = p.adjust(res$P)
  return(res)

}

correlations = lapply(c("SBS", "DBS", "ID"), score_versus_num_sigs)
do.call(rbind, correlations) -> correlations

correlations %>%
  ungroup() %>%
  huxtable::as_huxtable() %>%
  set_align("center") %>%
  set_wrap(TRUE) %>%
  insert_row(
    paste("Estimated coefficient and P values by linear regression",
          "with the R lm() function (Combined_Score ~ num_signatures).",
          "FDR is the",
          "Benjamini-Hochberg false discovery rate within",
          "each mutation type (SBS, DBS ID)."),
    fill = ""
  ) %>%
  insert_row(
    paste("Table S20: For each signature attribution approach,",
          " Combined Score as a function of",
          "the number of mutational signatures considered for",
          "attribution in a cancer type."),
    fill = ""
  ) %>%
  set_header_rows(c(1, 3), TRUE) %>%
  style_header_rows(bold = TRUE) %>%
  merge_across(1:2, 1:ncol(.)) %>%
  set_row_height(row = 1:2, value = 1.8 / nrow(.)) %>%
  set_row_height(row = 3, value = 1.4 / nrow(.)) %>%
  set_valign("middle") %>%
  quick_xlsx(file = "output_for_paper/table_s20_score_versus_number_of_sigs.xlsx")

