library(ggplot2)
library(ggbeeswarm)
library(gridExtra)
library(dplyr)

do_plot_ranks_by_tool = FALSE

mean_ranks_by_tool = function(mutation_type) {
  ranktable = fread(
    paste0("output_for_paper/tools_ranked_in_each_cancer_type_",
    mutation_type,
    ".csv"))

  
  ranktable = group_by(ranktable, Tool)
  
  means <- ranktable %>%
    group_by(Tool) %>%
    summarize(mean_rank = mean(`Rank within cancer type`)) %>%
    ungroup() %>%
    arrange(mean_rank)
  
  means %>%
    mutate(`Cancer type` = "", rank = "") %>%
    relocate(`Cancer type`, Tool, rank, mean_rank) %>%
    rename(`Mean rank` = mean_rank) %>%
    huxtable::as_huxtable() %>%
    huxtable::set_align("centre") %>%
    huxtable::set_number_format(NA) %>%
    quick_xlsx(file = 
                 file.path("output_for_paper",
                           paste0("add_this_manually_to_sup_table_for_",
                                  mutation_type, "_rank.xlsx")))
  
  if (do_plot_ranks_by_tool) {
    
    ggplot(ranktable, aes(x = Tool, y = `Rank within cancer type`, group = Tool)) +
      geom_beeswarm(aes(color = `Cancer type`), size = 2) +
      theme_minimal() +
      labs(
        title = paste("Ranks for", mutation_type, "within each cancer type"),
        x = "Tool",
        y = "Rank",
        color = "Cancer Type"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "grey90", size = 0.5)
      ) +
      scale_y_reverse(breaks = 1:13)
  }
}

if (do_plot_ranks_by_tool ) {
  # This plot is not used in the paper
  grid.arrange(
    plot_ranks_by_tool("SBS"),
    plot_ranks_by_tool("DBS"),
    plot_ranks_by_tool("ID"),
    ncol = 1, 
    nrow = 3  )
}
