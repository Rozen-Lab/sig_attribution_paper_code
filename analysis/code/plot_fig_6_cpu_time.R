# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())

source("analysis/code/common_utils.R")
library(ggplot2)
library(gridExtra)

cpu_barplot <- function(cpu_seconds_file, 
                        mutation_type,  # Temporarily broken
                        title = "") {
  cpu_time <- data.table::fread(cpu_seconds_file)
  # browser()
  cpu_time = 
    dplyr::mutate(cpu_time, 
                  log_cpu_seconds = 
                    log10(cpu_seconds + 1))
  colnames(cpu_time)[1] <- "Tool"
  
  tools_to_plot = tools_to_plot_and_order(mutation_type)
  
  cpu_time = dplyr::filter(cpu_time, Tool %in% tools_to_plot) # Need to deal with changed tools names
  
  # cpu_time <- change_tool_names(cpu_time)
  cpu_time$Tool = pretty_tool_names(cpu_time$Tool)
  
  tools = tools_to_plot_and_order("SBS") 
  custom_colors1 = 
    rev(RColorBrewer::brewer.pal(12, "Set3"))
  custom_colors = c(custom_colors1, "#999933")
  names(custom_colors) = pretty_tool_names(tools)
  
  # browser()
  # Reorder the tools by cpu.hloours
  cpu_time$Tool <- factor(cpu_time$Tool, levels = cpu_time$Tool[order(cpu_time$log_cpu_seconds)])
  
  # Create the bar plot
  ggplot(cpu_time, aes(x = Tool, y = log_cpu_seconds, fill = Tool)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = custom_colors) +
    # scale_y_log10(limits = c(0, 450)) +  # Log10 scale for y-axis
    labs(x = "Tool", y = "Log base 10 (CPU seconds + 1)", title = title) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
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
    ylim(0, 6.5) -> plot_object
  
  return(plot_object)
}

sbs <-
  cpu_barplot(file.path(global_summary_data_root,
                        "total_cpu_seconds_SBS.csv"),
              mutation_type = "SBS",
              title = "SBS"
  )


dbs <-
  cpu_barplot(file.path(global_summary_data_root, 
                        "total_cpu_seconds_DBS.csv"),
              mutation_type = "DBS",
              title = "DBS"
  )


id <-
  cpu_barplot(file.path(global_summary_data_root, 
                        "total_cpu_seconds_ID.csv"),
              mutation_type = "ID",
              title = "ID"
  )

ggplot2::ggsave(
  filename = file.path(global_output_for_paper, "fig_6_cpu_time.pdf"),
  plot = gridExtra::marrangeGrob(grobs = list(sbs, dbs, id),
                                 nrow = 1, ncol = 3),
  width = 8.2677, height = 5, units = "in"
)

