library(data.table)
source("analysis/code/analysis_utils.R")
source("common_code/common_utils.R")


# Used in analysis for sup table 5
get_missed_sig_summaries <-
  function(gt_exp_file, inferred_exp_files, tools, cancer_type) {
    retval <- lapply(seq_along(tools), FUN = function(index) {
      tmp <- get_missed_sig_summary(
        gt_exp_file = gt_exp_file,
        inferred_exp_file = inferred_exp_files[index],
        tool = tools[index],
        cancer_type = cancer_type
      )
      return(tmp)
    })
    
    df <- retval[[1]]
    
    for (i in 2:length(retval)) {
      df <- merge(x = df, y = retval[[i]])
    }
    sig_order <- gtools::mixedsort(df$gt_sig_id)
    df2 <- df[match(sig_order, df$gt_sig_id), ]
    
    gt_sig_info <-
      get_gt_sig_count(gt_exp_file = gt_exp_file, cancer_type = cancer_type)
    gt_sig_info2 <- gt_sig_info[df2$gt_sig_id, ]
    
    summary_table <- merge(x = gt_sig_info2, y = df2)
    summary_table2 <-
      summary_table[match(sig_order, summary_table$gt_sig_id), ]
    return(summary_table2)
  }




get_missed_sig_summary <-
  function(gt_exp_file, inferred_exp_file, tool, cancer_type) {
    gt_exp <- mSigTools::read_exposure(gt_exp_file)
    samples_one_type <-
      grep(pattern = cancer_type, x = colnames(gt_exp), value = TRUE)
    gt_exp_one_type <- gt_exp[, samples_one_type]
    inferred_exp <- mSigTools::read_exposure(inferred_exp_file)
    inferred_exp_one_type <- inferred_exp[, samples_one_type]
    
    missed_gt_sigs <-
      lapply(samples_one_type, FUN = function(sample_name) {
        gt_exp <-
          mSigAct:::RemoveZeroActivitySig(gt_exp_one_type[, sample_name,
                                                          drop = FALSE
          ])
        inferred_exp <-
          mSigAct:::RemoveZeroActivitySig(inferred_exp_one_type[, sample_name,
                                                                drop = FALSE
          ])
        missed_sigs <- setdiff(rownames(gt_exp), rownames(inferred_exp))
        gt_sig_exp <- gt_exp[missed_sigs, ]
        df <- data.frame(
          gt_sig_id = missed_sigs,
          gt_sig_exp = gt_sig_exp
        )
        rownames(df) <- NULL
        return(df)
      })
    
    missed_sig_table <- do.call(rbind, missed_gt_sigs)
    missed_sig_summary <- missed_sig_table %>%
      dplyr::group_by(gt_sig_id) %>%
      dplyr::summarise(
        missed_times = n()
      )
    colnames(missed_sig_summary)[2] <-
      paste0(tool, "_", colnames(missed_sig_summary)[2])
    return(missed_sig_summary)
  }






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
  "analysis/summary/SBS/assessment_each_sample_SBS.csv"
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

writexl::write_xlsx(
  x = list(`Table S5` = summary_table),
  path = file.path(global_output_for_paper, 
                   "sup_table_s5_syn_SBS_skin_missed_sig.xlsx"),
  col_names = TRUE, format_headers = TRUE)

write_csv(summary_table, 
          file.path(global_output_for_paper, 
                    "sup_table_s5_syn_SBS_skin_missed_sig.csv"),)