# Run this script with the top level directory as the working directory
stopifnot(basename(getwd()) == "sig_attribution_paper_code")
rm(list = ls())
source("analysis/code/generic_gather_and_print_stats.R")
source("analysis/code/gather_and_save_cpu_time.R")

base_sup_table = 6

generic_gather_and_print_stats("SBS", sup_table_number = base_sup_table) # tables S6, s7, s8, s9, s10

# Table s11 is generated by sourcing skin_missed_sig_sup_table.R

generic_gather_and_print_stats("DBS", sup_table_number = base_sup_table + 6) # tables s12, s13, s14, s15

generic_gather_and_print_stats("ID", sup_table_number = base_sup_table + 10) # tables s16, s17, s18, s19

sup_table_number = base_sup_table + 14 # tables s20, s21, s22
for (mutation_type in c("SBS", "DBS", "ID")) {
  gather_and_save_cpu_time(mutation_type, sup_table_number)
  sup_table_number = sup_table_number + 1
}

# Table S23 and sup figs S7, S8, S9
# are generated by sourcing cpu_versus_num_sigs.R

# Be sure to run mean_ranks_by_tool.R and put the mean ranks 
# in each rank sup table (sup tables s8, s13, s16 for SBS, DBS, ID, respectively)
