library(data.table)
source("analysis/code/analysis_utils.R")

sbs_gt_file <- "synthetic_data/SBS/ground.truth.syn.exposures.csv"
inferred_exp_files <-
  list.files(
    path = "analysis/raw_output/SBS/",
    pattern = "inferred_exposures.csv",
    full.names = TRUE,
    recursive = TRUE
  )
tools <- basename(sub("/syn.*", "", inferred_exp_files))

assesment_by_sample_file <-
  "analysis/summary/SBS/syn/assessment_each_sample.csv"
assesment_by_sample <- data.table::fread(assesment_by_sample_file)
tool_order <- get_tool_order(assesment_by_sample)
index_order <- match(tool_order, tools)

inferred_exp_files2 <- inferred_exp_files[index_order]
tools2 <- tools[index_order]
cancer_type <- "Skin-Melanoma"

summary_table <-
  get_missed_sig_summaries(
    gt_exp_file = sbs_gt_file,
    inferred_exp_files = inferred_exp_files2,
    tools = tools2,
    cancer_type = cancer_type
  )
colnames(summary_table)[1:2] <-
  c("Signature", "Number of tumors with signature")
colnames(summary_table) <-
  gsub(
    pattern = "_missed_times", replacement = "",
    x = colnames(summary_table)
  )


output_home <- "output_for_paper/"
writexl::write_xlsx(
  x = list(`Table S5` = summary_table),
  path = file.path(output_home, "sup_table_s5_syn_SBS_skin_missed_sig.xlsx"),
  col_names = TRUE, format_headers = TRUE
)
