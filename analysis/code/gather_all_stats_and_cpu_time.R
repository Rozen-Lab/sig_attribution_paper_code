# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/generic_gather_and_print_stats.R")
source("analysis/code/gather_and_save_cpu_time.R")

generic_gather_and_print_stats("SBS", sup_table_number = 5) # tables S5, s6, s7

generic_gather_and_print_stats("DBS", sup_table_number = 10) # tables s10, s11, s12

generic_gather_and_print_stats("ID", sup_table_number = 13) # tables s13, s14, s15

sup_table_number = 16 # tables s16, s17, s18
for (mutation_type in c("SBS", "DBS", "ID")) {
  gather_and_save_cpu_time(mutation_type, sup_table_number)
  sup_table_number = sup_table_number + 1
}

# Table S19 is generated from cpu_versus_num_sigs.R
