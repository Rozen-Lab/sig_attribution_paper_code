library(data.table)
source("analysis/code/analysis_utils.R")

sbs_file <- "analysis/summary/SBS/syn/assessment_each_sample.csv"
sbs_dt <- data.table::fread(sbs_file)
sbs_tool_order <- get_tool_order(sbs_dt)
sbs_test1 <- 
  wilcox.test(x = sbs_dt[Tool == "pasa",]$Combined, 
              y = sbs_dt[Tool == "sigpro",]$Combined)
sbs_p_value1 <- sbs_test1$p.value
sbs_test2 <- 
  wilcox.test(x = sbs_dt[Tool == "pasa",]$Combined, 
              y = sbs_dt[Tool == "fitms",]$Combined)
sbs_p_value2 <- sbs_test2$p.value

dbs_file <- "analysis/summary/DBS/syn/assessment_each_sample.csv"
dbs_dt <- data.table::fread(dbs_file)
dbs_tool_order <- get_tool_order(dbs_dt)
dbs_test1 <- 
  wilcox.test(x = dbs_dt[Tool == "pasa",]$Combined, 
              y = dbs_dt[Tool == "fitms",]$Combined)
dbs_p_value1 <- dbs_test1$p.value

id_file <- "analysis/summary/ID/syn/assessment_each_sample.csv"
id_dt <- data.table::fread(id_file)
id_tool_order <- get_tool_order(id_dt)
id_test1 <- 
  wilcox.test(x = id_dt[Tool == "pasa",]$Combined, 
              y = id_dt[Tool == "fitms",]$Combined)
id_p_value1 <- id_test1$p.value

