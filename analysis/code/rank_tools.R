library(tidyverse)
library(huxtable)
library(openxlsx)

rank_tools = function(mutation_type, sup_table_number, 
                      summary_stats, # Summary stats over all cancer types together
                      summary_stats_by_cancer_type)  {
  # browser()

  if (mutation_type == "SBS") {
    # Used to generate a supplementary table that ranks calls to fitms
    # called with different values of the "rare signature threshold.
    summary_stats %>% 
      mutate(V1 = NULL, base_tool = sub("_.*", "", Tool)) %>%
      filter(base_tool == "fitms") %>%
      separate(Tool, into = c("fitms", "threshold"), sep = "_", convert = TRUE) %>%
      arrange(base_tool, desc(m.Combined)) %>%
      mutate(rank = row_number(), base_tool = NULL) %>%
      relocate(rank, .after = threshold) %>% 
      select(2:4) %>%
      arrange(desc(m.Combined)) %>%
      rename(`Mean (Combined Score)` = m.Combined,
             Rank = rank, 
             `Rare signature threshold` = threshold)  -> 
      fitms_rank_by_rare_threshold
  }

  summary_stats_by_cancer_type %>%
    filter(Tool %in% global_raw_tools_to_plot) %>%
    group_by(cancer.type) %>%
    arrange(cancer.type, desc(m.Combined)) %>%
    mutate(rank_in_cancer_type = min_rank(desc(m.Combined))) %>%  # Ensures tied ranks
    relocate(cancer.type) %>%
    relocate(rank_in_cancer_type, .after = Tool) %>%
    select(cancer.type, Tool, rank_in_cancer_type, m.Combined) %>%
    rename(`Cancer type` = cancer.type, 
           `Rank within cancer type` = rank_in_cancer_type) %>%
    mutate(Tool = pretty_tool_names(Tool),
           `mean (Combined Score)` = m.Combined,
           .keep = "unused") ->
    tools_ranked_in_each_cancer_type

  
  f2 = function(var) {
    origname = deparse(substitute(var))
    outfile = file.path(
      global_output_for_paper, 
      paste0(origname, "_", mutation_type, ".csv"))
    # browser()
    data.table::fwrite(var, outfile)
  }
  
  if (mutation_type == "SBS") {
    write_sup_table_as_excel(fitms_rank_by_rare_threshold, 
                             my_sup_table_number = as.numeric(sup_table_number) + 2, 
                             mutation_type = mutation_type,
                             fitms = TRUE)
  }
  write_sup_table_as_excel(tools_ranked_in_each_cancer_type, 
                           my_sup_table_number = sup_table_number, 
                           mutation_type = mutation_type)
  f2(tools_ranked_in_each_cancer_type)
}


write_sup_table_as_excel = function(var, my_sup_table_number, mutation_type, fitms = FALSE) {
  if (fitms) {
    headingtext = 
      paste0("Table S", my_sup_table_number, 
             ": FitMS mean Combined Score for ",
             "different thresholds for rare ",
             "signatures")
  } else {
    headingtext = 
      paste0("XXXXXTable S", my_sup_table_number, 
             ": Signature attribution tool rank within each cancer type for ",
             mutation_type, " signatures",
             " XXXXXX Important manually add contents of ",
             "add_this_manually_to_sup_table_for_",
             mutation_type, "_rank.xlsx")
  }
  
  origname = deparse(substitute(var))
  outfile = file.path(
    global_output_for_paper, 
    paste0("table_s", my_sup_table_number, "_", origname, "_", mutation_type, ".xlsx"))
  # browser()
  var %>%
    huxtable::as_huxtable() %>%
    set_align("center") %>%
    insert_row(headingtext, fill = "") %>%
    set_valign("middle") %>%
    merge_cells(1, 1:ncol(.)) %>%
    set_header_rows(1:2, TRUE) %>%
    style_headers(bold = TRUE) %>%
    set_row_height(1:2, 2 / nrow(.)) %>%
    set_number_format(NA) ->
    xvar
  
  xxvar = as_Workbook(xvar)
  numstyle = createStyle(numFmt = "0.000")
  numeric_cols <- which(sapply(var, is.numeric))  # Identify numeric columns
  addStyle(xxvar, sheet = 1,
           style = numstyle,
           rows = 2:nrow(xvar),
           cols = 4,
           gridExpand = TRUE)
  setColWidths(xxvar, 1, 1:4, widths = c(25, 15, 15, 25))
  saveWorkbook(xxvar, file = outfile, overwrite = TRUE)
  
}


rank_tools_multiple_measures = function(summary_stats) {
  summary_stats %>%
    dplyr::filter(Tool %in% global_raw_tools_to_plot) %>%
    mutate(Tool = pretty_tool_names(Tool)) -> nexttable
  cols_to_analyze = 
    c("m.Combined", "sd.Combined", "m.F1", "m.sum_sens_spec", "m.SMD")
  nexttable = nexttable[ , c("Tool", cols_to_analyze)]
  for (colname in cols_to_analyze) {
    # browser()
    newcolname = paste0("rank_by_", colname)
    nexttable[ , newcolname] = min_rank(desc(pull(nexttable, colname)))
    dplyr::relocate(nexttable, newcolname, .after = colname) -> nexttable
  }

  # browser()
  maxrank_sd = max(nexttable$rank_by_sd.Combined)
  maxrank_SMD = max(nexttable$rank_by_m.SMD)
  nexttable$rank_by_sd.Combined = 1 + maxrank_sd - nexttable$rank_by_sd.Combined
  nexttable$rank_by_m.SMD       = 1 + maxrank_SMD - nexttable$rank_by_m.SMD
  colnames(nexttable) = c(
    "Tool",
    "Mean Combined Score",
    "Rank by mean Combined Score",
    "SD of Combined Score",
    "Rank by SD of Combined Score",
    "Mean F1",
    "Rank by mean F1",
    "Mean recall + specificity",
    "Rank by mean recall + specificity",
    "Mean scaled Manhattan distance",
    "Rank by mean scaled Manhattan distance"
  )
  return(nexttable) 
}



