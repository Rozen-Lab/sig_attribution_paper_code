gather_cpu_info = function(mutation_type) {

  output_home <- file.path("analysis/raw_output", mutation_type)
  total_cores <- parallel::detectCores()
  cores_to_use <- total_cores / 2
    
  list1  <-
    get_total_cpu_seconds(
      ouput_dir = output_home,
      mc_cores = cores_to_use
    )
  
  total_cpu_seconds          = list1$total_cpu_seconds
  cpu_seconds_by_cancer_type = list1$cpu_seconds_by_cancer_type

  if (FALSE) {
  elapsed_time <-
    get_elapsed_time(
      ouput_dir = output_home,
      mc_cores = cores_to_use
    )
  }
  
  # Get the tool order
  output_dir_syn <- file.path("analysis/summary", mutation_type, "syn")
  summary_stats <-
    read.csv(file = file.path(output_dir_syn, "all_summary_stats.csv"))
  summary_stats2 <- dplyr::arrange(summary_stats, desc(med.Combined))
  tool_order <- summary_stats2$Tool
  
  if (FALSE) { # Skipping odering rows for now 2024 09 10
    to2 = 1:length(tool_order)
    
    browser()
    
    names(to2) = tool_order
    df = dplyr::arrange(total_cpu_seconds, to2[tool_name])
  }
  

  write_1_table = function(value, basename) {
    fpath = file.path(output_dir_syn, basename)
    message("Writing ", fpath)
    write.csv(
      x = value,
      file = fpath,
      row.names = FALSE
    )
  }
  
  write_1_table(total_cpu_seconds, "total_cpu_seconds.csv")
  write_1_table(cpu_seconds_by_cancer_type, "cpu_seconds_by_cancer_type.csv")
  
  if (FALSE) {
    elapsed_time2 <- elapsed_time[tool_order]
    df2 <- data.frame(
      tool = names(elapsed_time2),
      elapsed_time_seconds = elapsed_time2
    )
    write.csv(
      x = df2,
      file = file.path(output_dir_syn, "elapsed_time_seconds.csv"),
      row.names = FALSE
    )
  }
}


get_total_cpu_seconds <- function(ouput_dir, mc_cores = 1) {
  rds_files <-
    list.files(
      path = ouput_dir,
      full.names = TRUE, recursive = TRUE,
      pattern = "^time_used.Rds"
    )
  
  # browser()
  total_cpu_seconds <- 
    lapply(rds_files, 
           function(file) {
             tool_name <- basename(sub("/syn.*", "", file))
             timing <- readRDS(file)
             return(
               tibble(tool_name = tool_name, total_cpu_sec = sum(timing[c(1, 2, 4, 5)])))
           })
  total_cpu_seconds = do.call(rbind, total_cpu_seconds)    
             
  rds_by_cancer_type = 
    list.files(
      path = ouput_dir,
      full.names = TRUE, recursive = TRUE,
      pattern = "^time_by_cancer_type.Rds"
    )
  
  process_cpu_by_cancer_type = function(in_rds_path) {
    tool_name <- basename(sub("/syn.*", "", in_rds_path))
    by_cancer_type_cpu <-readRDS(in_rds_path)
    xx = lapply(
      names(by_cancer_type_cpu), 
      function(cancer_type) {
        return(tibble(
          tool_name   = tool_name,
          cancer_type = cancer_type, 
          cpu_seconds = sum((by_cancer_type_cpu[[cancer_type]])[c(1, 2, 4, 5)])))
      })
    xx = do.call(rbind, xx)
    # browser()
    return(xx)
  }
  
  cpu_by_cancer_type = 
    do.call(
      rbind,
      lapply(rds_by_cancer_type, process_cpu_by_cancer_type))
  
  top_level_dirs <- list.dirs(path = ouput_dir, recursive = FALSE)
  msa_dirs <- grep(pattern = "msa", x = top_level_dirs, value = TRUE)
  
  if (length(msa_dirs) > 0) {
    msa_cpu_hours <- sapply(msa_dirs, FUN = function(msa_dir) {
      get_msa_total_cpu_hours(
        msa_output_dir = msa_dir,
        mc_cores = mc_cores
      )
    })
    msa_tool_names <- basename(sub("/syn.*", "", msa_dirs))
    names(msa_cpu_hours) <- msa_tool_names
    msa_cpu_seconds <- msa_cpu_hours * 60 * 60
    # browser()
    msa_table = tibble(tool_name = names(msa_cpu_seconds),
                       total_cpu_sec = msa_cpu_hours)
    total_cpu_seconds <- rbind(total_cpu_seconds, msa_table)
  } else {
    message("No MSA results")
  }
  dplyr::group_by(cpu_by_cancer_type, tool_name) %>% 
    dplyr::summarise(sum_of_cpu_by_cancer_type = sum(cpu_seconds)) ->
    sum_cpu_by_cancer_type
  
  # browser()
  new_total_cpu = dplyr::full_join(total_cpu_seconds, sum_cpu_by_cancer_type)

  return(list(total_cpu_seconds          = new_total_cpu, 
              cpu_seconds_by_cancer_type = cpu_by_cancer_type))
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


