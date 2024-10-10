# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
library(data.table)
library(tidyverse)
output_for_paper_dir = "output_for_paper"

## Results section for SBS

sbs_summary = fread(
  file.path(output_for_paper_dir,
            "SBS/all_summary_stats_SBS.csv"))
arrange(sbs_summary, desc(m.Combined)) %>% head() %>% select(Tool, m.Combined)

sbs_file <- file.path(output_for_paper_dir, 
                      "SBS/assessment_each_sample_SBS.csv")
sbs_dt <- data.table::fread(sbs_file)
sbs_test1 <- 
  wilcox.test(x = sbs_dt[Tool == "pasa",]$Combined, 
              y = sbs_dt[Tool == "musical",]$Combined)
sbs_p_value1 <- sbs_test1$p.value
sbs_test2 <- 
  wilcox.test(x = sbs_dt[Tool == "pasa",]$Combined, 
              y = sbs_dt[Tool == "fitms_0.010",]$Combined)
sbs_p_value2 <- sbs_test2$p.value

dplyr::filter(sbs_summary, Tool %in% c("pasa", "musical", "fitms_0.010")) %>%
         dplyr::select(Tool, m.Combined)

## Results section for DBS
dbs_summary = fread(
  file.path(output_for_paper_dir,
            "DBS/all_summary_stats_DBS.csv"))
arrange(dbs_summary, desc(m.Combined)) %>% head() %>% select(Tool, m.Combined)

dbs_file <- file.path(output_for_paper_dir, 
                      "DBS/assessment_each_sample_DBS.csv")
dbs_dt <- fread(dbs_file)
dbs_test1 <- 
  wilcox.test(x = dbs_dt[Tool == "pasa",]$Combined, 
              y = dbs_dt[Tool == "fitms_0.010",]$Combined)
dbs_p_value1 <- dbs_test1$p.value
dbs_test2 <- 
  wilcox.test(x = dbs_dt[Tool == "fitms_0.010",]$Combined, 
              y = dbs_dt[Tool == "musical",]$Combined)
dbs_p_value2 <- dbs_test2$p.value

## Manuscript Results section for ID

id_summary = fread(
  file.path(output_for_paper_dir,
            "ID/all_summary_stats_ID.csv")) %>%
  arrange(desc(m.Combined)) %>% head(20) %>% 
  select(Tool, m.Combined) %>% print()


id_file <- file.path(output_for_paper_dir, 
                     "ID/assessment_each_sample_ID.csv")
id_dt <- data.table::fread(id_file)

wilcox.test(x = id_dt[Tool == "pasa",]$Combined, 
            y = id_dt[Tool == "fitms_0.060",]$Combined) %>% `$`(p.value)

wilcox.test(x = id_dt[Tool == "fitms_0.060",]$Combined, 
            y = id_dt[Tool == "musical",]$Combined) %>% `$`(p.value)
