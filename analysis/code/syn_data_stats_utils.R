library(tibble)
library(openxlsx)
source("analysis/code/get_all_input.R")
source("analysis/code/common_utils.R")

get_stats <- function(spectra, cancer_type) {
  mut_counts <- colSums(spectra)
  df <- tibble::tibble(
    "Cancer type" = cancer_type,
    "Mean" = mean(mut_counts),
    "Median" = median(mut_counts),
    "SD" = sd(mut_counts)
  )
  return(df)
}

get_stats_by_mut_type <- function(mut_type) {
  syn_data <- get_all_input(mut_type)
  all_type_stats <-
    get_stats(syn_data$all_spectra,
      cancer_type = "All cancer types"
    )
  by_type_spectra <- syn_data$spectra_list
  by_type_stats <-
    lapply(names(by_type_spectra), FUN = function(cancer_type) {
      one_type_spectra <- by_type_spectra[[cancer_type]]
      one_type_stats <- get_stats(
        spectra = one_type_spectra,
        cancer_type = cancer_type
      )
    })
  by_type_stats2 <- do.call(rbind, by_type_stats)
  all_stats <- rbind(all_type_stats, by_type_stats2)
  return(all_stats)
}

get_stats_all_syn_data <- function() {
  mut_types <- c("SBS", "DBS", "ID")

  all_syn_data_stats <- lapply(mut_types, FUN = function(mut_type) {
    return(get_stats_by_mut_type(mut_type))
  })

  names(all_syn_data_stats) <- mut_types
  wb <- createWorkbook()

  addWorksheet(wb, "table_s3")

  title_style <-
    createStyle(
      fontSize = 11, textDecoration = "bold",
      halign = "left"
    )
  header_style <-
    createStyle(
      fontSize = 11, textDecoration = "bold",
      halign = "center", valign = "center"
    )
  content_style <-
    createStyle(
      fontSize = 11, halign = "center",
      valign = "center", numFmt = "0.00"
    )

  writeData(wb, "table_s3", all_syn_data_stats[["SBS"]],
    startRow = 3, startCol = 1
  )
  writeData(wb, "table_s3", all_syn_data_stats[["DBS"]][, 2:4],
    startRow = 3, startCol = 5
  )
  writeData(wb, "table_s3", all_syn_data_stats[["ID"]][, 2:4],
    startRow = 3, startCol = 8
  )
  writeData(wb, "table_s3", "SBS", startRow = 2, startCol = 2, colNames = FALSE)
  writeData(wb, "table_s3", "DBS", startRow = 2, startCol = 5, colNames = FALSE)
  writeData(wb, "table_s3", "ID", startRow = 2, startCol = 8, colNames = FALSE)

  mergeCells(wb, "table_s3", cols = 2:4, rows = 2)
  mergeCells(wb, "table_s3", cols = 5:7, rows = 2)
  mergeCells(wb, "table_s3", cols = 8:9, rows = 2)

  writeData(wb, "table_s3", "Table S3: Statistics for synthetic data",
    startRow = 1, startCol = 1, colNames = FALSE
  )
  addStyle(wb, "table_s3", title_style,
    rows = 1,
    cols = 1, gridExpand = TRUE
  )
  addStyle(wb, "table_s3", header_style,
    rows = 2:3,
    cols = 1:10, gridExpand = TRUE
  )
  addStyle(wb, "table_s3", content_style,
    rows = 4:13,
    cols = 1:10, gridExpand = TRUE
  )
  setColWidths(wb, "table_s3", cols = 1, widths = 20)
  saveWorkbook(wb, file.path(global_output_for_paper,
    "table_s3_synthetic_data_stats.xlsx"),
    overwrite = TRUE
  )
}
