library(tidyverse)
library(huxtable)
library(openxlsx)

rank_tools = function(mutation_type, sup_table_number, summary_stats, summary_stats_by_cancer_type)  {
  # browser()
  if (FALSE) {
    df = data.table::fread(
      file.path(global_output_for_paper,
                paste0("summary_stats_by_cancer_type_",
                       mutation_type, ".csv")))

  xx = data.table::fread(
    file.path(global_output_for_paper, 
              paste0("all_summary_stats_",
                     mutation_type, ".csv")))
  } else {
    df = summary_stats_by_cancer_type
    xx = summary_stats
  }

  if (mutation_type == "SBS") {
     # Used to generate the file "best_each_tool_<mutation_type>.csv"
    xx %>% 
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

  df %>%
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

  
  f1 = function(var, my_sup_table_number, fitms = FALSE) {
    if (fitms) {
      headingtext = 
        paste0("Table S", my_sup_table_number, 
               ": FitMS mean Combined Score for ",
               "different thresholds for rare ",
               "signatures")
    } else {
      headingtext = 
        paste0("Table S", my_sup_table_number, 
               ": Signature attribution tool rank within each cancer type for ",
               mutation_type, " signatures")
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
  
  f2 = function(var) {
    origname = deparse(substitute(var))
    outfile = file.path(
      global_output_for_paper, 
      paste0(origname, "_", mutation_type, ".csv"))
    # browser()
    data.table::fwrite(var, outfile)
  }
  
  if (mutation_type == "SBS") {
    f1(fitms_rank_by_rare_threshold, my_sup_table_number = "8", fitms = TRUE)
  }
  f1(tools_ranked_in_each_cancer_type, my_sup_table_number = sup_table_number)
  f2(tools_ranked_in_each_cancer_type)
}

