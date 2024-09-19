# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())

library(rempsyc)
library(flextable)

create_comparison_table = function(mutation_type, table_nunber) {
  browser()
  file <-
    file.path("analysis/summary", mutation_type, "syn/assessment_each_sample.csv")
  dt <- data.table::fread(file)
  dt$one_minus_SMD <- 1 - dt$SMD
  
  median_combined <-
    aggregate(dt$Combined,
              list(dt$Tool),
              FUN = median
    )
  tools <-
    median_combined[order(median_combined$x, decreasing = TRUE), ]$Group.1
  metrics <-
    c("Combined", "one_minus_SMD", "prec", "sens")
  
  format_number <- function(x) {
    format(round(x, 2), nsmall = 2)
  }
  
  med_iqr_info <- lapply(tools, FUN = function(tool) {
    statistics <- character()
    for (metric in metrics) {
      tmp <- summary(dt[Tool == tool, get(metric)])
      info <-
        paste0(
          format_number(tmp["Median"]), " (",
          format_number(tmp["1st Qu."]),
          " to ", format_number(tmp["3rd Qu."]), ")"
        )
      statistics <- c(statistics, info)
    }
    df <- data.frame(statistics)
    colnames(df) <- tool
    rownames(df) <- metrics
    return(df)
  })
  
  med_iqr_df <- do.call(cbind, med_iqr_info)
  
  metric_df <-
    data.frame(Metric = c(
      "Combined", "1 - SMD", 
      "Precision", "Recall"
    ))
  med_iqr_df <- cbind(metric_df, med_iqr_df)
  my_table <- rempsyc::nice_table(
    data = med_iqr_df,
    title = c(paste("Table", table_number), 
              paste("Performance of attribution tool on synthetic",
              mutation_type, "data")),
    note = c(
      "Each cell in the table indicates the median with lower quartile and upper quartile in the bracket",
      "SMD: Scaled Manhattan Distance"
    )
  )
  
  output_home <- file_path("analysis/summary", mutation_type, "syn/ID")
  table_name <- paste0("comparison_table_", mutation_type, ".docx")
  flextable::save_as_docx(my_table,
                          path = file.path(output_home, table_name)
  )
}

create_comparison_table("ID", "3")
