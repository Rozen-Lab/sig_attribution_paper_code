# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())

source("analysis/code/plot_functions.R")

## CPU PLOT #####################

cpu_barplot <- function(cpu_seconds_file, main = "",
                        ylim = c(0, 450),
                        mutation_type) { # Temporarily broken
  cpu_time <- data.table::fread(cpu_seconds_file)
  
  cpu_time = 
    dplyr::mutate(cpu_time, 
                  total_cpu_hours = 
                    (sum_of_cpu_by_cancer_type / 60) / 60)
  colnames(cpu_time)[1] <- "Tool"
  
  cpu_time = dplyr::filter(cpu_time, Tool %in% raw_tools_to_plot(mutation_type)) # Need to deal with changed tools names
  
  cpu_time <- change_tool_names(cpu_time)
  
  if (any(duplicated(cpu_time$Tool))) {
    browser()
  }
  cpu_time$Tool <- factor(cpu_time$Tool, levels = cpu_time$Tool)
  
  
  
  if (main == "") {
    title <- "Total CPU Time of Different Tool"
  } else {
    title <- main
  }
  
  plot_object <-
    ggplot(cpu_time, aes(x = Tool, y = total_cpu_hours, fill = Tool)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = get_colors(mutation_type)) +
    labs(
      title = title,
      x = "Tool",
      y = "Total CPU Time (hours)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      legend.position = "none",
      plot.margin = margin(
        t = 10,
        r = 15,
        b = 5,
        l = 30
      )
    ) +
    ylim(ylim[1], ylim[2])
  return(plot_object)
}



sbs_cpu_seconds_file <-
  "analysis/summary/SBS/total_cpu_seconds_SBS.csv"
sbs_cpu_plot_object <-
  cpu_barplot(sbs_cpu_seconds_file,
              main = "Synthetic SBS",
              ylim = c(0, 450),
              mutation_type = "SBS"
              
  )

dbs_cpu_seconds_file <-
  "analysis/summary/DBS/total_cpu_seconds_DBS.csv"
dbs_cpu_plot_object <-
  cpu_barplot(dbs_cpu_seconds_file,
              main = "Synthetic DBS",
              ylim = c(0, 450),
              mutation_type = "DBS"
  )

id_cpu_seconds_file <-
  "analysis/summary/ID/total_cpu_seconds_ID.csv"
# browser()
id_cpu_plot_object <-
  cpu_barplot(id_cpu_seconds_file,
              main = "Synthetic ID",
              ylim = c(0, 450),
              mutation_type = "ID"
  )

output_home <- "output_for_paper/"

ggplot_to_pdf(
  plot_objects = c(
    list(sbs_cpu_plot_object),
    list(dbs_cpu_plot_object),
    list(id_cpu_plot_object)
  ),
  file = "figure_6_cpu_time_three_types.pdf",
  nrow = 3, ncol = 3,
  width = 8.2677, height = 11.6929, units = "in"
)

