library(ggplot2)
dd = data.table::fread("analysis/summary/SBS/assessment_each_sample_SBS.csv")
check_LLH = function(dd) {

  dd_summary <- dd[, .(mean_val = mean(multiLLH), median_val = median(multiLLH)), by = Tool]
  
  
  dd_truncated <- dd[multiLLH >= -5000]
  
  ggplot(dd_truncated, aes(x = multiLLH, fill = Tool)) +
    geom_histogram(binwidth = 0.5, color = "black", alpha = 0.7) +  # Adjust binwidth as needed
    facet_wrap(~ Tool, scales = "free_y",ncol = 1) +  # Separate histograms by Tool
    theme_minimal() +
    labs(
      title = "Histogram of multiLLH Values for Each Tool (Truncated at -3000)",
      x = "multiLLH",
      y = "Frequency"
    ) +
    geom_vline(data = dd_summary, aes(xintercept = mean_val), color = "blue", linetype = "solid", size = 1) +  # Solid blue line for mean
    geom_vline(data = dd_summary, aes(xintercept = median_val), color = "red", linetype = "dashed", size = 1) +  # Dashed red line for median
    theme(
      legend.position = "none"  # Remove legend if you don't need it
    )
  
  }

check_LLH(dd[dd$Tool %in% c("sigpro", "mp")])
