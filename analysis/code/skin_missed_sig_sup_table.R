library(data.table)
library(mSigTools)
library(huxtable)
library(openxlsx)
library(tidyverse)
# source("analysis/code/analysis_utils.R")
source("analysis/code/common_utils.R")


get_gt_sig_count <- function(gt_exp_table, cancer_type) {
  samples_one_type <-
    grep(pattern = cancer_type, x = colnames(gt_exp_table), value = TRUE)
  gt_exp_one_type <-
    mSigAct:::RemoveZeroActivitySig(gt_exp_table[, samples_one_type])
  tmp <- apply(X = gt_exp_one_type, MARGIN = 1, FUN = function(counts) {
    return(length(counts[counts > 0]))
  })
  df <- data.frame(gt_sig_id = names(tmp), gt_sig_count = tmp)
  return(df)
}


# Use in analysis for a sup table
get_missed_sig_summaries_SBS <-
  function(cancer_type, only_some_tools = TRUE) {
    gt_exp_table = mSigTools::read_exposure(
      "synthetic_data/SBS/ground.truth.syn.exposures.csv")
    inferred_exp_table = fread(
      "output_for_paper/all_inferred_exposures_SBS.csv")
    tools <- tools_to_plot_and_order("SBS")
    if (only_some_tools) {
     # tools = setdiff(tools, c("pasa", "sigfit"))
      tools = setdiff(tools, "sigfit")
    }
    # browser()
        
    retval <- lapply(tools, FUN = function(tool) {
      tmp <- get_missed_sig_summary_one_tool(
        gt_exp_table = gt_exp_table,
        inferred_exp_table = inferred_exp_table,
        tool = tool,
        cancer_type = cancer_type
      )
      return(tmp)
    })

    outdf <-
      get_gt_sig_count(gt_exp_table = gt_exp_table, cancer_type = cancer_type)
    
    for (retdf in retval) {
      outdf = full_join(outdf, retdf)
    }
    # browser()   
    outdf2 = mutate(outdf, across(everything(), ~replace_na(., 0)))

    colnames(outdf2)[1:2] <-
      c("Signature", "Number of tumors with signature")
    
    return(outdf2)
  }


get_missed_sig_summary_one_tool <-
  function(gt_exp_table, inferred_exp_table, tool, cancer_type) {

    samples_one_type <-
      grep(pattern = cancer_type, x = colnames(gt_exp_table), value = TRUE)
    gt_exp_one_type <- gt_exp_table[, samples_one_type]

    dplyr::filter(inferred_exp_table, 
                  Tool == tool, Sample.ID %in% samples_one_type) %>%
      mutate(Sample.ID = NULL, Tool = NULL) %>% t ->
      inferred_exp_one_type
    colnames(inferred_exp_one_type) = samples_one_type
 
    
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
        if (tool == "sigfit" && length(missed_sigs) > 0) message("sigfit missed a signature")
        gt_sig_exp <- gt_exp[missed_sigs, ]
        df <- data.frame(
          gt_sig_id = missed_sigs,
          gt_sig_exp = gt_sig_exp
        )
        rownames(df) <- NULL

        return(df)
      })
    # if (tool == "sigfit") browser()
    missed_sig_table <- do.call(rbind, missed_gt_sigs)
    missed_sig_summary <- missed_sig_table %>%
      dplyr::group_by(gt_sig_id) %>%
      dplyr::summarise(
        missed_times = n()
      )
    colnames(missed_sig_summary)[2] <- pretty_tool_names(tool)
    return(missed_sig_summary)
  }


sum_table <- get_missed_sig_summaries_SBS("Skin-Melanoma")
# summary_table = get_missed_sig_summaries_SBS("Skin-Melanoma", FALSE)


generate_huxtable = function(mytable) {
  huxtable::as_huxtable(mytable) %>% 
    insert_row("",
               "Number of synthetic tumors with signature",
               "Number of tumors in which the signature was missed", 
               fill = "",
               after = 0) %>%
    merge_cells(1, 3:ncol(.)) %>% 
    insert_row(
      paste("For each tool, red indicates the most-often-missed signature,",
            "and orange indicates the second-most-often missed signature.",
            "sigfit was omitted because it assigns every possible signature to every",
            "Skin-Melanoma tumor."), 
      fill = "", after = 0) %>%
    merge_cells(1, 1:ncol(mytable)) %>%
    insert_row("Table S11: In Skin-Melanoma, the most commonly missed signatures were SBS1 and SBS5",
               fill = "", after = 0) %>% ## Change table number below as well
    merge_cells(1, 1:ncol(mytable)) %>%
    set_header_rows(c(1, 3,4), TRUE) %>%
    style_headers(bold = TRUE) %>%
    set_bold(c(1, 3:nrow(.)), 1) %>%
    set_align(3:nrow(.), 1:ncol(.), "center") %>%
    set_bottom_border(c(2, 4, 15),  everywhere, value = 1) %>%
    set_bottom_border(3, 3:ncol(.), value = 1) %>%
    merge_cells(3:4, 2) %>%
    set_valign("middle") -> tmp_table
  
  for (index in 3:ncol(mytable)) {
    message(index)
    rows = which(mytable[ , index] == max(mytable[ , index], na.rm = TRUE))
    tmp_table = set_background_color(tmp_table, rows + 4, index, "red1")
    
  }
  
  for (index in 3:ncol(mytable)) {
    rows = get_second_max_rows(mytable[ , index])
    tmp_table = set_background_color(tmp_table, rows + 4, index, "tan1")
    
  }
  
  return(tmp_table)
}




# Function to get the row(s) with the second largest value for each column
get_second_max_rows <- function(column) {
  # Remove duplicates and get the second largest value
  unique_values <- sort(unique(column), decreasing = TRUE)
  if (length(unique_values) < 2) {
    return(NA) # If there's no second largest value
  }
  second_largest <- unique_values[2]
  # Find the rows where the column equals the second largest value
  which(column == second_largest)
}


table1 = generate_huxtable(sum_table)

quick_xlsx( # Change table number above on title line of spreadsheet as well
  table1,
  file = 
    file.path(global_output_for_paper, 
              "table_s11_skin_missed_sig.xlsx"))


if (FALSE) {
  TT = as_tibble(sum_table)
  sum_table2 = TT %>% mutate(across(3:ncol(.), ~ .x / TT[[2]]))
  
  table2 = generate_huxtable(sum_table2)
  quick_xlsx(table2, 
             file = file.path(global_output_for_paper, "do_not_use_for_now.xlsx"))
}

# Check assertion that for Skin-Melanoma sigfit estimated that every signature appeared in every tumor

sigfit <- read_csv("analysis/raw_output/SBS/sigfit/syn/inferred_exposures.csv")
sigs = sigfit$...1
sigfit = data.frame(sigfit[ , grepl("Skin", colnames(sigfit))])
dim(sigfit)
rownames(sigfit) = sigs
sigfit = sigfit[sum_table[, 1], ] # only the signatures in Skin-Melanoma
any(sigfit < 0.5)
