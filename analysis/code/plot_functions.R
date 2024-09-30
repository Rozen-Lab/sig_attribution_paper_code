
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

source("analysis/code/common_utils.R")

## GLOBAL VARIABLES AND COMMON FUNCTIONS #####################

old_old_custom_colors <-
  c(
    "pasa" = "#E69F00", "sigpro" = "#009E73", "fitms" = "#F0E442",
    "mp" = "#0072B2", "msa" = "#CC79A7", "msa_opt" = "#56B4E9",
    "msa_thresh_x10" = "#999933", "msa_thresh_x100" = "#D55E00",
    "msa_thresh_x1000" = "#999999"
  )

old_custom_colors <-
  c(
    "PASA" = "#E69F00", "SigPro" = "#009E73", "FitMS" = "#F0E442",
    "MP" = "#0072B2", "MSA" = "#CC79A7", "MSA_opt" = "#56B4E9"
  )


if (length(global_raw_tools_to_plot) > 12) {
  extra = length(global_raw_tools_to_plot) - 12
  custom_colors1 = 
    rev(RColorBrewer::brewer.pal(12, "Set3"))
  custom_colors = c(custom_colors1, custom_colors1[1:extra])
} else {
  custom_colors =
    rev(RColorBrewer::brewer.pal(length(global_raw_tools_to_plot), "Set3"))
}

names(custom_colors) = global_raw_tools_to_plot

ggplot_to_pdf <-
  function(plot_objects, file, nrow, ncol, width, height, units, out_directory = plot_output_directory) {
    ggplot2::ggsave(
      filename = file.path(out_directory, file),
      plot = gridExtra::marrangeGrob(plot_objects,
                                     nrow = nrow, ncol = ncol,
                                     layout_matrix = matrix(seq_len(nrow * ncol),
                                                            nrow = nrow, ncol = ncol,
                                                            byrow = TRUE
                                     )
      ),
      width = width, height = height, units = units
    )
  }

scale_fn <- function(x) {
  sprintf("%.2f", x)
}


## PLOTS WITH ALL CANCER TYPES COMBINED #####################

plot_all_cancer_types_merged = function(mut_type, fig_num) {
  
  data_home <- file.path("analysis/summary", mut_type)
  orig_indata <-
    data.table::fread(
      file.path(data_home, 
                paste0("assessment_each_sample_", mut_type, ".csv")))
  indata <- orig_indata # change_tool_names(orig_indata)
  
  tools_to_plot = global_raw_tools_to_plot

  plot_objects <-
    boxplots_combined_cancer_types(
      assessment_by_sample = indata,
      tools_to_plot = tools_to_plot,
      mut_type = mut_type
    )
  
  ggplot_to_pdf(
    plot_objects = plot_objects,
    file = paste0("figure_", fig_num, "_", mut_type, "_all_cancer_types.pdf"),
    nrow = 4, ncol = 2,
    width = 8.2677, height = 11.6929, units = "in"
  )
}


boxplots_combined_cancer_types <-
  function(assessment_by_sample, main = "", 
           tools_to_plot = NULL,
           mut_type) {
    # assessment_by_sample is a big data.table
    # with one row for each run of one tool on 
    # one input spectrum and the columns being
    # various measures such MD (Manbhattan distance),
    # sens, prec, Combined, etc.

    if (!is_null(tools_to_plot)) {
      assessment_by_sample = 
        dplyr::filter(assessment_by_sample, Tool %in% tools_to_plot)
    }
    # browser()
    tool_order <- get_tool_order(assessment_by_sample) 
    # tool_order is a vector of character strings
    
    assessment_by_sample$cancer.type <- "All cancer types"
    # browser()
    df <-
      data.frame(dplyr::mutate(assessment_by_sample,
                               one_minus_smd = 1 - SMD))

    df$Tool <- factor(df$Tool, levels = tool_order)
    # browser()
    measures <- c("Combined", "one_minus_smd", "prec", "sens", "spec", 
                  "scaled_L2", "KL", "multiLLH")
    ylabs <- global_ylabs[measures]
    
    my_plots <- 
      lapply(seq_along(measures),
             FUN = function(index) {
               # Make a plot for each measure (e.g. "Combined", "Precision", etc.)
               one_boxplot_combined_cancer_types(
                 measure = measures[index],
                 df = df,
                 xlab = "Tool", 
                 legend_position = "none",
                 ylab = ylabs[index], main = main,
                 mut_type = mut_type
               )
             })
    
    return(my_plots)
  }


one_boxplot_combined_cancer_types <-
  function(measure, df, xlab, legend_position,
           ylab = measure, main = NULL, mut_type) {
    # browser()
    measure <- sym(measure)
    xlab <- sym(xlab)
    if (is.null(main)) {
      main = mut_type
    }
    plot_object <-
      ggplot(df, aes(x = !!xlab, y = !!measure, fill = Tool)) +
      geom_boxplot(outlier.size = 0.1, outlier.color = "grey") +
      stat_summary(
        fun = mean, geom = "point", shape = 18,
        size = 2, color = "red", fill = "red"
      ) +
      scale_fill_manual(values = custom_colors) +
      scale_x_discrete(labels = function(x) sapply(x, function(toolname) (global_tools_to_plot[toolname]))) +
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
      ggtitle(mut_type) +
      ylab(ylab)
    # message("main is ", main)
    
    if (measure == "KL") { 
      return(plot_object+ scale_y_reverse())
    } 
    if (measure == "multiLLH") {
      if (mut_type == "SBS") {
        return(plot_object + coord_cartesian(ylim = c(-1500, -200)))
      } else if (mut_type == "DBS") {
        return(plot_object + coord_cartesian(ylim = c(-500, 0)))
      } else if (mut_type == "ID") {
        return(plot_object + coord_cartesian(ylim = c(-750, 0)))
      } else {
        stop("Unknown mutation type: ", mut_type)
      }
    } 
    return(plot_object  + coord_cartesian(ylim = c(0, NA)))
    
  }


## PLOTS BY EACH CANCER TYPE #####################

boxplots_by_cancer_type <-
  function(assessment_by_sample, 
           main = "", 
           tools_to_plot = NULL,
           measures =
             c("Combined", "one_minus_smd", "prec", "sens"),
           mut_type
  ) {
    # Boxplots of many measures by individual cancer type
    
    ylabs = global_ylabs[measures]
    
    
    # assessment_by_sample is a big data.table
    # with one row for each run of one tool on 
    # one input spectrum and the columns being
    # various measures such MD (Manhattan distance),
    # sens, prec, Combined, etc.
    # browser()
    
    # browser()
    if (!is_null(tools_to_plot)) {
      assessment_by_sample = 
        dplyr::filter(assessment_by_sample, Tool %in% tools_to_plot)
    }
    
    
    tool_order <- get_tool_order(assessment_by_sample)
    
    df <-
      data.frame(dplyr::mutate(assessment_by_sample,
                               one_minus_smd = 1 - SMD
      ))
    
    df$Tool <- factor(df$Tool, levels = tool_order)
    
    
    by_type_plots <- list()
    
    for (index in seq_along(measures)) {
      retval <-
        one_boxplot_by_cancer_type(
          measure = measures[index],
          df = df,
          xlab = "cancer.type", legend_position = "none",
          ylab = ylabs[index], main = main,
          last_plot = (index == length(measures)),
          mutation_type = mut_type
        )
      by_type_plots <- c(by_type_plots, list(retval))
    }
    
    return(by_type_plots)
  }

plot_by_cancer_type_sup_fig = function(mutation_type, sup_fig_num) {
  
  
  data_home <- file.path("analysis/summary", mutation_type)
  indata <-
    data.table::fread(
      file.path(data_home, 
                paste0("assessment_each_sample_", mutation_type, ".csv")))
  
  # indata <- change_tool_names(indata)
  # browser()
  tools_to_plot = global_raw_tools_to_plot
  if (mutation_type != "DBS") {
    tools_to_plot = tools_to_plot[-global_msa_index]
  }
  
  plot_some = function(measures, basename) {
    # browser()
    plot_objects1 <-
      boxplots_by_cancer_type(
        assessment_by_sample = indata,
        tools_to_plot = tools_to_plot,
        measures = measures,
        mut_type = mutation_type
      )
    
    ggplot_to_pdf(
      plot_objects = plot_objects1,
      file = basename,
      nrow = 2, ncol = 1,
      width = 11.6929, 
      height = 8.2677,
      units = "in"
    )
  }
  

  filepat = paste0("sup_fig_", sup_fig_num, "YYY_", mutation_type, "_by_cancer_type.pdf")
  
  plot_some(c("Combined", "one_minus_smd"),
            basename = gsub("YYY", "AB", filepat))
  
  plot_some(c("prec", "sens"),
            basename = gsub("YYY", "CD", filepat))
  
  plot_some(c("spec", "scaled_L2"),
            basename = gsub("YYY", "EF", filepat))
  
  plot_some(c("KL", "multiLLH"),
            basename = gsub("YYY", "GH", filepat))
}

one_boxplot_by_cancer_type <-
  function(measure, df, xlab, legend_position,
           ylab = measure, mutation_type = NULL,
           main = NULL, last_plot = FALSE) {
    # browser()
    last_plot = TRUE
    if (!is.null(mutation_type)) {
      main = mutation_type
    }
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
      scale_fill_manual(values = custom_colors) +
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
      ylab(ylab) # +
      # scale_y_continuous(labels = scale_fn) + coord_cartesian(ylim = c(0, NA))
    
    if (last_plot) {
      plot_object <-
        plot_object +
        theme(
          axis.text.x = element_text(angle = 0, vjust = 0.5),
          legend.position = "bottom",
          # axis.text.y = element_blank(),
          # axis.title.y = element_blank(),
          plot.margin = margin(
            t = 0,
            r = 30,
            b = 50,
            l = 30
          )
        ) + 
        guides(fill = guide_legend(nrow = 1))
    }
    
    if (measure == "KL") { 
      return(plot_object+ scale_y_reverse())
    } 
    if (measure == "multiLLH") {
      if (mutation_type == "SBS") {
        return(plot_object + coord_cartesian(ylim = c(-1500, -200)))
      } else if (mutation_type == "DBS") {
        return(plot_object + coord_cartesian(ylim = c(-500, 0)))
      } else if (mutation_type == "ID") {
        return(plot_object + coord_cartesian(ylim = c(-750, 0)))
      } else {
        stop("Unknown mutation type: ", mut_type)
      }
    } 
    return(plot_object  + coord_cartesian(ylim = c(0, NA)))
  }



## CPU PLOT #####################

cpu_barplot <- function(cpu_seconds_file, main = "", ylim = c(0, 450)) {
  cpu_time <- data.table::fread(cpu_seconds_file)
  
  cpu_time = 
    dplyr::mutate(cpu_time, total_cpu_hours = (sum_of_cpu_by_cancer_type / 60) / 60)
  colnames(cpu_time)[1] <- "Tool"
  cpu_time <- change_tool_names(cpu_time)
  cpu_time = dplyr::filter(cpu_time, Tool %in% global_tools_to_plot)
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
    scale_fill_manual(values = custom_colors) +
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

