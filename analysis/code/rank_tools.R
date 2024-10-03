library(tidyverse)

rank_tools = function(mut_type, view = FALSE)  {
  # browser()
  df = data.table::fread(
    file.path("analysis/summary/",
              mut_type, 
              paste0("summary_stats_by_cancer_type_",
                     mut_type, ".csv")))
  
  xx = data.table::fread(
    file.path("analysis/summary/",
              mut_type, 
              paste0("all_summary_stats_",
                     mut_type, ".csv")))
  
  xx %>% mutate(V1 = NULL, base_tool = sub("_.*", "", Tool)) %>%
    relocate(base_tool, Tool) %>%
    group_by(base_tool) %>% 
    arrange(base_tool, desc(m.Combined)) %>%
    mutate(tool_rank = row_number()) %>%
    relocate(tool_rank, .after = Tool) %>% 
    filter(tool_rank == 1) %>%
    arrange(desc(m.Combined)) -> 
    best_each_tool
  
  select(best_each_tool, Tool) %>%
    left_join(df) %>% 
    ungroup() %>%
    mutate(V1 = NULL, base_tool = NULL) %>% 
    group_by(cancer.type) %>%
    arrange(cancer.type, desc(m.Combined)) %>%
    mutate(rank_in_cancer_type = row_number()) %>%
    relocate(cancer.type) %>%
    relocate(rank_in_cancer_type, .after = Tool) %>%
    select(cancer.type, Tool, rank_in_cancer_type, m.Combined) ->
    best_tools_ranked_in_each_cancer_type
  
  filter(best_tools_ranked_in_each_cancer_type, rank_in_cancer_type == 1) ->
    best_tool_each_cancer_type
  
  filter(best_tools_ranked_in_each_cancer_type, Tool == "pasa") ->
    pasa_rank_each_cancer_type
  
  f1 = function(var) {
    origname = deparse(substitute(var))
    if (view) {
      View(var, paste0(mut_type, "_", origname))
    }
    print(origname)
    outfile = file = file.path(
      "analysis/summary", mut_type, 
      paste0(origname, "_", mut_type, ".csv"))
    # browser()
    data.table::fwrite(var, outfile)
  }
  
  f1(best_each_tool)
  f1(best_tools_ranked_in_each_cancer_type)
  f1(best_tool_each_cancer_type)
  f1(pasa_rank_each_cancer_type)
}
