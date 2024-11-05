library(dplyr)
library(rvest)
library(xml2) # read_html for MSA output parsing
library(tidyverse)
library(huxtable)

gather_and_save_cpu_time = function(mutation_type, sup_table_number) {
  
  raw_cpu_data_dir <- file.path("analysis/raw_output", mutation_type)
  total_cores <- parallel::detectCores()
  cores_to_use <- total_cores / 2
  if (Sys.info()["sysname"] == "Windows") cores_to_use = 1
  # cores_to_use = 1 # For debugging
  
  rds_by_cancer_type = 
    list.files(
      path = raw_cpu_data_dir,
      full.names = TRUE, recursive = TRUE,
      pattern = "^time_by_cancer_type.Rds"
    )
  
  process_cpu_by_cancer_type = function(in_rds_path) {
    tool_name <- basename(sub("/syn.*", "", in_rds_path))
    by_cancer_type_cpu <-readRDS(in_rds_path)
    
    xx = lapply(
      names(by_cancer_type_cpu), 
      function(cancer_type) {
        if (FALSE && tool_name == "sigspack") {
          browser()
        }
        return(dplyr::tibble(
          tool_name   = tool_name,
          cancer_type = cancer_type, 
          cpu_seconds = 
            sum((by_cancer_type_cpu[[cancer_type]])[c(1, 2, 4, 5)],
                na.rm = TRUE)))
      })
    xx = do.call(rbind, xx)
    # browser()
    return(xx)
  }
  
  cpu_by_cancer_type = 
    do.call(
      rbind,
      lapply(rds_by_cancer_type, process_cpu_by_cancer_type))
  top_level_dirs <- list.dirs(path = raw_cpu_data_dir, recursive = FALSE)
  msa_dirs <- grep(pattern = "msa", x = top_level_dirs, value = TRUE)
  
  if (length(msa_dirs) > 0) {
    # Ignore unpruned output
    msa_dirs <- grep(pattern = "unpruned", x = msa_dirs, 
                     invert = TRUE, value = TRUE)
    # Only look at default thresholds
    msa_dirs = grep("default", msa_dirs, value = TRUE)
    
    msa_cpu_time <- lapply(msa_dirs, FUN = function(msa_dir) {
      get_msa_cpu_time(
        msa_output_dir = msa_dir,
        mc_cores = cores_to_use
      )
    })

    msa_cpu_sec_by_type <- do.call(rbind, lapply(msa_cpu_time, `[[`, 2))
    
    cpu_by_cancer_type <- rbind(cpu_by_cancer_type, msa_cpu_sec_by_type)
  } else {
    message("No MSA results")
  }
  
  cpu_by_cancer_type  %>%
    dplyr::group_by(tool_name) %>% 
    dplyr::summarise(
      cancer_type = "Total all cancer types",
      cpu_seconds = sum(cpu_seconds)) ->
    total_all_cancer_types
  
  bind_rows(cpu_by_cancer_type, total_all_cancer_types) %>%
    arrange(tool_name, cancer_type) %>%
    filter(tool_name %in% global_raw_tools_to_plot) ->
    cpu_by_cancer_type
  
  write_1_table = function(value, basename) {
    fpath = file.path("output_for_paper", basename)
    message("Writing ", fpath)
    write.csv(
      x = value,
      file = fpath,
      row.names = FALSE
    )
  }
  
  # We need this file as input for plotting
  write_1_table(total_all_cancer_types, 
                paste0("total_cpu_seconds_", mutation_type, ".csv"))

  # We need this file as input for test association between 
  # number of signatures and time complexity of attribution
  write_1_table(cpu_by_cancer_type, 
                paste0("table_s",
                       sup_table_number, 
                       "_cpu_by_cancer_type_", mutation_type, ".csv"))
  
  sum_columns = 1 + which(grepl("Total", cpu_by_cancer_type$cancer_type, fixed = TRUE))
  cpu_by_cancer_type %>%
    mutate(tool_name = pretty_tool_names(tool_name)) %>%
    rename(Tool = tool_name, 
           `Cancer type` = cancer_type,
           `CPU seconds` = cpu_seconds) %>%
    huxtable::as_huxtable() %>%
    set_align("center") %>%
    set_valign("middle") %>%
    huxtable::set_bold(
      row = sum_columns,
      col = 1:ncol(.),
      TRUE
    ) %>%
    insert_row(
      paste0("Table S", sup_table_number, ": ",
             "CPU seconds for each cancer type for ",
             mutation_type, " signartures"),
      fill = "",
      after = 0
    ) %>% 
    merge_cells(1, 1:ncol(.)) %>%
    set_header_rows(1:2, TRUE) %>%
    style_header_rows(bold = TRUE) %>%
    set_row_height(1:2, 2 / nrow(.)) %>%
    quick_xlsx(
      file = file.path(
        "output_for_paper",
        paste0("table_s",
               sup_table_number, 
               "_cpu_by_cancer_type_", mutation_type, ".xlsx")))
}


get_msa_cpu_hours <- function(html_file) {
  html <- read_html(html_file)
  
  all_text <- html %>%
    html_nodes("div") %>%
    html_text()
  
  tmp <- grep(pattern = "CPU-Hours", x = all_text, value = TRUE)[1]
  tmp2 <- unlist(stringr::str_split(string = tmp, pattern = "\n"))
  index <- grep(pattern = "CPU-Hours", x = tmp2)
  cpu_hours <- as.numeric(tmp2[index + 1])
  return(cpu_hours)
}

get_msa_cpu_time <- function(msa_output_dir, mc_cores) {
  tool_name <- basename(msa_output_dir)
  msa_time_files <- list.files(
    path = msa_output_dir,
    full.names = TRUE, recursive = TRUE,
    pattern = "^MSA-nf_report.html"
  )
  msa_time_files2 <-
    grep(
      pattern = "/raw/", x = msa_time_files, invert = TRUE,
      value = TRUE
    )
  cancer_types <- basename(dirname(dirname(msa_time_files2)))
  msa_cpu_hours <-
    parallel::mclapply(msa_time_files2,
                       FUN = get_msa_cpu_hours,
                       mc.cores = mc_cores
    )
  msa_cpu_hours_raw <- 
    dplyr::tibble(cpu_hours =  unlist(msa_cpu_hours),
                  cancer_type = cancer_types)
  msa_cpu_hours_by_type <- msa_cpu_hours_raw %>%
    dplyr::group_by(cancer_type) %>%
    dplyr::summarise(cpu_hours_by_type = sum(cpu_hours))
  
  total_cpu_sec <- 
    dplyr::tibble(tool_name = tool_name, 
                  total_cpu_sec = sum(unlist(msa_cpu_hours)) * 60 * 60)
  cpu_sec_by_type <-
    dplyr::tibble(tool_name = tool_name, 
                  cancer_type = msa_cpu_hours_by_type$cancer_type,
                  cpu_seconds = msa_cpu_hours_by_type$cpu_hours_by_type * 60 * 60)
  return(list(total_cpu_sec = total_cpu_sec,
              cpu_sec_by_type = cpu_sec_by_type))
}

get_msa_total_cpu_hours <- function(msa_output_dir, mc_cores) {
  msa_time_files <- list.files(
    path = msa_output_dir,
    full.names = TRUE, recursive = TRUE,
    pattern = "^MSA-nf_report.html"
  )
  msa_time_files2 <-
    grep(
      pattern = "/raw/", x = msa_time_files, invert = TRUE,
      value = TRUE
    )
  
  msa_cpu_hours <-
    parallel::mclapply(msa_time_files2,
                       FUN = get_msa_cpu_hours,
                       mc.cores = mc_cores
    )
  
  total_cpu_hours <- sum(unlist(msa_cpu_hours))
  return(total_cpu_hours)
}

get_msa_total_elapsed_time <- function(msa_output_dir, mc_cores) {
  msa_time_files <- list.files(
    path = msa_output_dir,
    full.names = TRUE, recursive = TRUE,
    pattern = "^MSA-nf_report.html"
  )
  msa_time_files2 <-
    grep(
      pattern = "/raw/", x = msa_time_files, invert = TRUE,
      value = TRUE
    )
  
  msa_elapsed_time <-
    parallel::mclapply(msa_time_files2,
                       FUN = get_msa_elapsed_time,
                       mc.cores = mc_cores
    )
  
  msa_total_elapsed_time <- sum(unlist(msa_elapsed_time))
  return(msa_total_elapsed_time)
}

get_msa_elapsed_time <- function(html_file) {
  html <- read_html(html_file)
  
  all_text <- html %>%
    html_nodes("div") %>%
    html_text()
  
  tmp <- grep(pattern = "duration", x = all_text, value = TRUE)[1]
  tmp2 <- unlist(stringr::str_split(string = tmp, pattern = "\n"))
  tmp3 <- grep(pattern = "duration", x = tmp2, value = TRUE)
  minute_used <-
    as.numeric(na.omit(stringr::str_extract(
      string = tmp3,
      pattern = "(\\d)+(?=m)"
    )))
  seconds_used <-
    as.numeric(na.omit(stringr::str_extract(
      string = tmp3,
      pattern = "(\\d)+(?=s)"
    )))
  elapsed_seconds <- minute_used * 60 + seconds_used
  return(elapsed_seconds)
}


