
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(grid)
library(gridExtra)

source("analysis/code/common_utils.R")

## GLOBAL VARIABLES AND COMMON FUNCTIONS #####################

get_colors = function(mutation_type) {
  tools = tools_to_plot_and_order("SBS") 
  custom_colors1 = 
    rev(RColorBrewer::brewer.pal(12, "Set3"))
  custom_colors = c(custom_colors1, "#999933")
  names(custom_colors) = tools
  return(custom_colors)
}

ggplot_to_pdf <-
  function(plot_objects, file, nrow, ncol, width, height, units, 
           out_directory = plot_output_directory) {
    ggplot2::ggsave(
      filename = file.path(out_directory, file),
      plot = gridExtra::marrangeGrob(plot_objects,
                                     nrow = nrow, ncol = ncol,
                                     layout_matrix = matrix(seq_len(nrow * ncol),
                                                            nrow = nrow, ncol = ncol,
                                                            byrow = TRUE
                                     )
      ),
      width = width, height = height, units = units, onefile = TRUE
    )
  }

scale_fn <- function(x) {
  sprintf("%.2f", x)
}


## PLOTS WITH ALL CANCER TYPES COMBINED #####################

plot_all_cancer_types_merged = function(mutation_type, fig_num, sup_fig_num) {

  indata <-
    data.table::fread(
      file.path(global_summary_data_root,  
                paste0("stats_each_sample_", mutation_type, ".csv")))

  tools_to_plot = tools_to_plot_and_order(mutation_type) # raw_tools_to_plot(mutation_type)

  plot_objects <-
    boxplots_combined_cancer_types(
      assessment_by_sample  = indata,
      ordered_tools_to_plot = tools_to_plot,
      mutation_type         = mutation_type
    )
  
  ggplot_to_pdf(
    plot_objects = plot_objects[1:4],
    file = paste0("fig_", fig_num, "_", mutation_type, "_all_cancer_types.pdf"),
    nrow = 4, ncol = 2,
    width = 8.2677, height = 11.6929, units = "in"
  )
  
  ggplot_to_pdf(
    plot_objects = plot_objects[5:length(plot_objects)],
    file = paste0("sup_fig_S", sup_fig_num, "_", mutation_type, "_all_cancer_types.pdf"),
    nrow = 4, ncol = 2,
    width = 8.2677, height = 11.6929, units = "in"
  )
}


boxplots_combined_cancer_types <-
  function(assessment_by_sample, main = "", 
           ordered_tools_to_plot,
           mutation_type) {
    # assessment_by_sample is a big data.table
    # with one row for each run of one tool on 
    # one input spectrum and the columns being
    # various measures such MD (Manbhattan distance),
    # sens, prec, Combined, etc.
    
    assessment_by_sample = 
      dplyr::filter(assessment_by_sample, Tool %in% ordered_tools_to_plot)
    
    # browser()
    # tool_order <- get_tool_order(assessment_by_sample)
    # tool_order is a vector of character strings
    
    assessment_by_sample$cancer.type <- "All cancer types"
    # browser()
    df <-
      data.frame(dplyr::mutate(assessment_by_sample,
                               one_minus_smd = 1 - SMD))

    df$Tool <- factor(df$Tool, levels = ordered_tools_to_plot)
    # browser()
    # measures <- c("Combined", "one_minus_smd", "prec", "sens", "spec", 
    #              "scaled_L2", "KL")
    ylabs <- global_ylabs[global_measures]
    
    my_plots <- 
      lapply(seq_along(global_measures),
             FUN = function(index) {
               # Make a plot for each measure (e.g. "Combined", "Precision", etc.)
               one_boxplot_combined_cancer_types(
                 measure = global_measures[index],
                 df = df,
                 xlab = "Tool", 
                 legend_position = "none",
                 ylab = ylabs[index], 
                 main = main,
                 mutation_type = mutation_type
               )
             })
    
    return(my_plots)
  }


one_boxplot_combined_cancer_types <-
  function(measure, df, xlab, legend_position,
           ylab = measure, main = "", mutation_type) {

    measure <- sym(measure)
    xlab <- sym(xlab)
    plot_object <-
      ggplot(df, aes(x = !!xlab, y = !!measure, fill = Tool)) +
      geom_boxplot(outlier.size = 0.1, outlier.color = "grey") +
      stat_summary(
        fun = mean, geom = "point", shape = 18,
        size = 2, color = "red", fill = "red"
      ) +
      scale_fill_manual(values = get_colors(mutation_type)) +
      scale_x_discrete(
        labels = function(x)  { 
          sapply(x,
                 function(toolname) { 
                   return(pretty_tool_names(toolname))
                 })
          }
        ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = legend_position,
        plot.margin = margin(
          t = 5,
          r = 30,
          b = 5,
          l = 30
        )
      ) + 
      ggtitle(main) +
      ylab(ylab)
    
    if (measure == "KL") { 
      return(plot_object+ scale_y_reverse())
    } 
    if (measure == "Combined") {
      if (mutation_type %in% c("SBS", "ID")) {
        return(plot_object + coord_cartesian(ylim = c(1, 3)))
      } else {
        return(plot_object + coord_cartesian(ylim = c(0.5, 3)))
      }
    }
    return(plot_object  + coord_cartesian(ylim = c(0, NA)))
    
  }


## PLOTS BY EACH CANCER TYPE #####################


plot_by_cancer_type_sup_fig = function(mutation_type, sup_fig_num) {
  
  indata <-
    data.table::fread(
      file.path(global_summary_data_root,
                paste0("stats_each_sample_", mutation_type, ".csv")))
  
  tools_to_plot = tools_to_plot_and_order(mutation_type)

  assessment_by_sample = dplyr::filter(indata, Tool %in% tools_to_plot)
  
  pdf(file.path(
    plot_output_directory,
    paste0("sup_fig_S", sup_fig_num, "_", mutation_type, ".pdf")),
    width = 11.6929, 
    height = 8.2677)

  for (ii in 1:length(global_measures)) {
    
    grob = one_boxplot_by_cancer_type(
      measure = global_measures[[ii]],
      df = assessment_by_sample,
      xlab = "cancer.type", legend_position = "none",
      last_plot = TRUE,
      mutation_type = mutation_type,
      main = paste0("\n\nSupplementary Figure S",
                    sup_fig_num,
                    LETTERS[ii], ", ", 
                    global_ylabs[global_measures[ii]], 
                    " by cancer type for ", mutation_type)
    )

    print(grob)

  }
  dev.off()

}

one_boxplot_by_cancer_type <-
  function(measure, df, 
           xlab, legend_position,
           # ylab = measure, 
           mutation_type,
           main, 
           last_plot = FALSE) {
    
    ylab = global_ylabs[measure]
    
    
    df <-
      data.frame(dplyr::mutate(df,
                               one_minus_smd = 1 - SMD
      ))
    
    tool_order <- tools_to_plot_and_order(mutation_type) # get_tool_order(df)

    df$Tool <- factor(df$Tool, levels = tool_order)

    last_plot = TRUE

    measure <- sym(measure)
    xlab <- sym(xlab)
    plot_object <-
      ggplot(df, aes(x = !!xlab, y = !!measure, fill = Tool)) +
      geom_boxplot(outlier.size = 0.1, outlier.color = "grey") +
      # This puts a red dot in the middle of each Tool within
      # each cancer type
      stat_summary(
        fun = mean, geom = "point", shape = 18,
        size = 2, color = "red", position = position_dodge(0.75)
      ) +
      scale_fill_manual(values = get_colors(mutation_type), labels = function(x) pretty_tool_names(x)) +
      # scale_fill_discrete(labels = function(x) pretty_tool_names(x)) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(
          t = 0,
          r = 30,
          b = 0,
          l = 30
        )
      ) +
      ggtitle(main) +
      ylab(ylab) 
    
    if (last_plot) {
      plot_object <-
        plot_object +
        theme(
          axis.text.x = element_text(angle = 0, vjust = 0.5),
          legend.position = "bottom",
          legend.text = element_text(angle = 90, vjust = 1),
          plot.margin = margin(
            t = 0,
            r = 30,
            b = 50,
            l = 30
          )
        ) + 
        guides(fill = guide_legend(nrow = 1, label.position = "bottom"))
    }
    # browser()
    if (measure == "KL") { 
      if (mutation_type %in% c("SBS", "DBS")) {
        return(
          plot_object = plot_object + coord_cartesian(ylim = c(3, 0)) +
            scale_y_reverse()
          )
      } else {
        return(
          plot_object + coord_cartesian(ylim = c(2.5, 0)) +
            scale_y_reverse()
        )
      }
      return(plot_object+ scale_y_reverse())
    } 
    if (measure == "Combined Score") {
      if (mutation_type %in% c("SBS", "ID")) {
        return(plot_object + coord_cartesian(ylim = c(1, 3)))
      } else {
        return(plot_object + coord_cartesian(ylim = c(0.5, 3)))
      }
    }
    return(plot_object  + coord_cartesian(ylim = c(0, NA)))
  }
