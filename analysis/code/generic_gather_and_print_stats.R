# This file writes the following files
#
# all_inferred_exposures_<mutation_type>.csv
# stats_each_sample_<mutation_type>.csv
# summary_stats_<mutation_type>.csv
# tools_ranked_in_cancer_type_<mutation_type>.csv
#
# table_s<X>_summary_stats_<mutation_type>.
# table_s<X+1>_summary_stats_by_cancer_type_<mutation_type>.xlsx
# table_s<X+3>_ranks_by_multiple_measures_<mutation_type>.xlsx
# 
# and calls rank_tools() to write
#
# table_s<X+2>_tools_ranked_in_each_cancer_type_<mutation_type>.xlsx
#
# and, for SBS only (from_rank_tools())
#
# table_s<X+4>_fitms_rank_by_rare_threshold_SBS.xlsx

library(parallel)
library(mSigTools)
library(philentropy)
library(tidyverse)
library(huxtable)
source("analysis/code/common_utils.R")
source("analysis/code/get_all_input.R")
source("analysis/code/rank_tools.R")


generic_gather_and_print_stats = function(mutation_type, sup_table_number) {
  
  inferred_exp_output_files <- list.files(
    path = file.path("analysis/raw_output", mutation_type),
    full.names = TRUE, recursive = TRUE,
    pattern = "^inferred_exposures.csv"
  )
  tool_names <- basename(sub("/syn.*", "", inferred_exp_output_files))
  
  all_inferred_exposures <- lapply(inferred_exp_output_files, FUN = function(file) {
    return(mSigTools::read_exposure(file))
  })
  names(all_inferred_exposures) <- tool_names
  
  exp_tool_dfs <-
    lapply(names(all_inferred_exposures), 
           FUN = function(name) {
             exp_tool_df <- add_tool_etc(
               exposure = all_inferred_exposures[[name]],
               tool = name
             )
             return(exp_tool_df)
           })

  exposure_all <- do.call(dplyr::bind_rows, exp_tool_dfs)
  exposure_all[is.na(exposure_all)] <- 0

  data.table::fwrite(
    exposure_all,
    file.path(global_output_for_paper, 
              paste0("all_inferred_exposures_", mutation_type, ".csv")))
  
  compute_and_write_stats(
    exposure_all = exposure_all,
    mutation_type = mutation_type,
    sup_table_number = sup_table_number
  )
}

summarize_measures = function(tb) {
  # tb is tibble with one row for each sample.
  # The columns include the basic statistics
  # for each sample.
  #
  # This function summarize them based on groups
  # that are already present in tb.

  summarize(
    tb,
    
    m.Combined    = mean(Combined),
    med.Combined  = median(Combined),
    sd.Combined   = sd(Combined),
    
    m.F1          = mean(F1),
    med.F1        = median(F1),
    sd.F1         = sd(F1),
    
    m.SMD         = mean(SMD),
    med.SMD       = median(SMD),
    sd.SMD        = sd(SMD),
    
    m.sens        = mean(sens),
    med.sens      = median(sens),
    sd.sens       = sd(sens),
    
    m.prec        = mean(prec),
    med.prec      = median(prec),
    sd.prec       = sd(prec),
    
    m.spec        = mean(spec),
    med.spec      = median(spec),
    
    m.scaled_L2   = mean(scaled_L2),
    med.scaled_L2 = median(scaled_L2),
    
    m.KL          = mean(KL),
    med.KL        = median(KL),
    
    m.sum_sens_spec   = mean(sum_sens_spec),
    med.sum_sens_spec = median(sum_sens_spec)
    
  ) -> retval
  
  return(retval)
}

compute_and_write_stats <- 
  function(exposure_all,
           mutation_type,
           sup_table_number) {

  message("writing statistics to directory ", global_output_for_paper)

  
  stats <- all_measures(exposure_all,
                        mutation_type)
  # stats is a very long list. 
  # Each element of stats will be one row in the matrix s2 on the next line
  s2 <- t(matrix(unlist(stats), nrow = length(stats[[1]]), byrow = FALSE))
  colnames(s2) <- c(
    "Sample.ID", "Tool", "MD", "SMD",
    "sens", "prec", "F1", "Combined", "spec", 
    "scaled_L2", "KL", "sum_sens_spec")
  
  as_tibble(s2) %>%
    mutate(
      SMD = as.numeric(SMD),
      sens = as.numeric(sens),
      prec = as.numeric(prec),
      F1 = as.numeric(F1),
      Combined = as.numeric(Combined),
      spec = as.numeric(spec),
      scaled_L2 = 1 - as.numeric(scaled_L2),
      KL        = log2(as.numeric(KL) + 1),
      sum_sens_spec = as.numeric(sum_sens_spec),
      cancer.type = sub("::.*", "", Sample.ID)
    ) -> assessment_each_sample

  assessment_each_sample %>%
    group_by(Tool) %>% 
    summarize_measures() ->
    summary.stats
  
  summary.stats2 <- dplyr::arrange(summary.stats, desc(m.Combined))
  
  write.csv( # We need the csv file for downstream processing
    summary.stats2,
    file = file.path(
      global_output_for_paper, 
      paste0("summary_stats_", mutation_type, ".csv")))

  summary.stats2 %>%
    filter(Tool %in% global_raw_tools_to_plot) %>%
    mutate(Tool = pretty_tool_names(Tool)) %>%
    `colnames<-`(
      unlist(lapply(colnames(.), 
                    rename_col))) %>%
    huxtable::as_huxtable() %>%
    set_align("center") %>%
    set_valign("middle") %>%
    insert_row(
      paste0("Table S", sup_table_number,
             ": Statistics across all cancer types for ",
             mutation_type, " attributions"),
      fill = ""
    ) %>%
    merge_cells(1, 1:ncol(.)) %>%
    set_row_height(2, 2/nrow(.)) %>%
    set_header_rows(1:2, TRUE) %>%
    style_header_rows(bold = TRUE) %>%
    quick_xlsx(
      file = file.path(
        global_output_for_paper, 
        paste0("table_s", 
               as.character(sup_table_number),
               "_summary_stats_", mutation_type, ".xlsx")))

  assessment_each_sample %>%
    group_by(Tool, cancer.type) %>%
    summarize_measures() ->
    summary.stats.by.cancer.type

  nncol = ncol(summary.stats.by.cancer.type)
  colwidths = rep(1, nncol)  / nncol 
  colwidths[1] = 2.5 * colwidths[1]
  colwidths = colwidths / sum(colwidths)
  summary.stats.by.cancer.type %>%
    filter(Tool %in% global_raw_tools_to_plot) ->
    summary2
  gindex = which((as.integer(as.factor(sort(summary2$cancer.type))) %% 2) == 0)
  
  summary2 %>%
    group_by(cancer.type) %>%
    arrange(cancer.type, desc(m.Combined)) %>%
    relocate(cancer.type) %>%
    mutate(Tool = pretty_tool_names(Tool)) %>%
    `colnames<-`(
      unlist(lapply(colnames(.), 
                    rename_col))) %>%
    huxtable::as_huxtable() %>%
    set_align("center") %>%
    set_valign("middle") %>%
    set_wrap(TRUE) %>%
    set_background_color(
      gindex + 1,
      col = 1:ncol(.), 
      value = "grey90") %>%
    insert_row(
      paste0("Table S", sup_table_number + 1, 
             ": Statistics by cancer type for ",
             mutation_type, " attributions"),
      fill = ""
    ) %>%
    set_header_rows(1:2, TRUE) %>%
    style_header_rows(bold = TRUE) %>%
    merge_across(1, 1:ncol(.)) %>%
    set_row_height(2, 2.5 / nrow(.)) %>%
    set_col_width(colwidths) %>%
    quick_xlsx(
      file = file.path(
        global_output_for_paper, 
        paste0("table_s", as.character(sup_table_number + 1), 
               "_summary_stats_by_cancer_type_", mutation_type, ".xlsx")))

  data.table::fwrite(
    assessment_each_sample,
    file = file.path(
      global_output_for_paper,
      paste0("stats_each_sample_", mutation_type, ".csv")))
  
  rank_tools(mutation_type, 
             sup_table_number = as.character(sup_table_number + 2),
             summary_stats = summary.stats2, 
             summary_stats_by_cancer_type = summary2)
  
  next_sup_table_number = sup_table_number + 3
  ranks = rank_tools_multiple_measures(summary.stats2)
  huxtable::as_huxtable(ranks) %>%
    set_align("center") %>%
    set_valign("middle") %>%
    insert_row(
      paste0(
        "Table S", next_sup_table_number,
        ": Tool ranks for " ,
        mutation_type, 
        " attributions for several measures"),
      fill = "") %>%
    merge_cells(1, 1:ncol(.)) %>%
    set_header_rows(1:2, TRUE) %>%
    style_headers(bold = TRUE) %>%
    set_row_height(2, 3.5 / nrow(.)) %>%
    quick_xlsx(
      file = file.path(
        global_output_for_paper, 
        paste0("table_s", 
               as.character(next_sup_table_number),
               "_ranks_by_multiple_measures_", mutation_type, ".xlsx")))

}

rename_col = function(tt) {
  if (tt == 'cancer.type') {
    return("Cancer type")
  }
  if (tt == "Tool") {
    return(tt)   
  }
  
  ttt = strsplit(tt, ".", fixed = TRUE)
  t1 = ttt[[1]][1]
  t2 = ttt[[1]][2]
  if(t1 == "m") {
    t1 = "Mean"
  } else if (t1 == "med") {
    t1 = "Median"
  } else if (t1 == "sd") {
    t1 = "SD"
  }

  if (t2 == "Combined") {
    t2 = "Combined Score"
  } else if (t2 == "SMD") {
    t2 = "scaled Manhattan distance"
  } else if (t2 == "sens") {
    t2 = "recall"
    # t2 = "sensitivity"
  } else if (t2 == "prec") {
    t2 = "precision"
  } else if (t2 == "scaled_L2") {
    t2 = "scaled L2 divergence"
  } else if (t2 == "spec") {
    t2 = "specificity"
  } else if (t2 == "sum_sens_spec") {
    # browser()
    t2 = "recall + specificity"
  }
  return(paste(t1,t2))
}


# Compute all measures for each row (sample) in table xx,
# which contains 1 exposure vector per row (sample).
all_measures <- function(xx, # A data.frame containing expsures. 
                         #Samples are in rows and signatures are in columns
                         mutation_type) { 
  
  total_cores <- parallel::detectCores()
  mc_cores <- total_cores / 2
  if (Sys.info()["sysname"] == "Windows") mc_cores = 1
  
  sig_universe_mask_file = 
    file.path("synthetic_data",
              mutation_type, # {SBS, DBS, ID}
              "sig_universe_mask.csv")
  sig_universe_mask = mSigTools::read_exposure(sig_universe_mask_file)
  
  ground_truth_exposures <- get_ground_truth_exposure(mutation_type)
  all_input = get_all_input(mutation_type)

  expected_spectra = mSigAct::ReconstructSpectrum(
    sigs          = all_input$ground_truth_sigs,
    exp           = ground_truth_exposures,
    use.sig.names = TRUE
  )
  
  one.row <- function(ii) {
    # ii is a row index in xx
    rr <- xx[ii, ]
    # rr is a 1-row tibble
    sid <- rr$Sample.ID
    tool <- rr$Tool
    stopifnot(!is.null(sid))
    stopifnot(!is.null(tool))
    
    if (!(sid %in% colnames(ground_truth_exposures))) {
      sid2 = gsub("MSI.H", "MSI-H", sid, fixed = TRUE)
      sid3 = sub(".", "-", sid2, fixed = TRUE)
      sid4 = sub("..", "::", sid3, fixed = TRUE)
      if (!(sid4 %in% colnames(ground_truth_exposures))) {
        browser()
      } else {
        sid = sid4
      }
    }

    gt0 <- ground_truth_exposures[ , sid]
    
    one_sig_universe_mask = t(sig_universe_mask[, sid, drop = FALSE] == 1)
    # gt0 <- dplyr::mutate(gt0, Sample.ID = NULL, Tool = NULL) # DELETE ME, no column Sample.ID
    me0 <- dplyr::mutate(rr, Sample.ID = NULL, Tool = NULL)
    stopifnot(sort(colnames(me0)) == sort(colnames(one_sig_universe_mask)))
    if (ncol(me0) > length(gt0)) {
      # This is because two signatures, SBS24 and SBS39 did not
      # have any mutations in any tumor
      stopifnot(setequal(
        setdiff(colnames(me0), names(gt0)), c("SBS24", "SBS39")))
      gt0["SBS24"] = 0
      gt0["SBS39"] = 0
    }
    # Align the signature order in gt0 and me0 to the
    # signature order from one_sig_universe_mask
    gt0 = gt0[colnames(one_sig_universe_mask)]
    me0 = me0[ , colnames(one_sig_universe_mask)]
    stopifnot(names(gt0) == colnames(me0))

    # Only look at signatures that were offered in the universe of signatures
    # for this sample.
    
    gt = gt0[one_sig_universe_mask]
    me = me0[ , one_sig_universe_mask]
    
    # Sanity check
    anti_gt = gt0[!one_sig_universe_mask]
    if(sum(anti_gt) != 0) {
      if (!(sid %in% c("Liver-HCC::S.36", "Liver-HCC::S.51"))) {
        message("Ground truth has signature not in the sig universe mask")
        message("Sample = ", sid)
        print(one_sig_universe_mask)
        print(gt0)
      }
    }
    
    me = unlist(me)
    gt = unlist(gt)

    if (all(me == 0)) {
      KLQ = me
    } else {
      KLQ = me / sum(me)
    }
    KL = philentropy::kullback_leibler_distance(
      P = gt / sum(gt), # The actual, ground truth distribution exposures
      Q = KLQ, # Our estimated exposures
      testNA = FALSE, unit = "log2", epsilon = 0.001)
    inferred_spectrum = mSigAct::ReconstructSpectrum(
      sigs          = all_input$ground_truth_sigs,
      exp           = as.matrix(me),
      use.sig.names = TRUE
    )
    
    if (FALSE) {
      multiLLH = mSigAct:::LLHSpectrumMultinom(
        spectrum = inferred_spectrum,
        expected.counts = expected_spectra[ , sid]
      )
    }
    
    L2 <- philentropy::distance(
      rbind(gt, me), method = "euclidean", mute.message = TRUE)
    
    manhattan <- sum(abs(gt - me))
    
    # if (any(me > 0 & me < 0.5)) { # Which tools create tiny exposures...?
    #  message(rr$Tool)
    # }
    
    called.pos = which(me >= 0.5) # This stems from the numerical optimization in 
                                 # sigspak, sigest, and presumably mutsig, in which
                                 # the optimization does not quite reach 0 from
                                 # above or below
    called.neg = which(me < 0.5)
    stopifnot(length(me) == (length(called.pos) + length(called.neg)))
    
    real.pos <- which(gt > 0)
    real.neg <- which(gt <= 0)
    
    TP <- length(intersect(called.pos, real.pos))
    
    P <- length(real.pos)
    
    sens <- TP / P

    TN <- length(intersect(called.neg, real.neg))
    N <- length(real.neg)
    
    if (N == 0) {
      spec = 1
    } else {
      spec = TN / N
    }
    
    FP <- length(setdiff(called.pos, real.pos))
    
    FN <- length(setdiff(called.neg, real.neg))
    # FN is correct because called.neg and real.neg are inflated by the
    # the same signatures not in the signature universe.
    
    prec <- TP / (TP + FP)
    if ((TP + FP) == 0) {
      prec <- 1
    }
    
    F1 <- (2 * TP) / (2 * TP + FP + FN)
    
    sumgt <- sum(gt)
    if (sumgt < 1) {
      browser()
    }
    scaled.manhattan <- manhattan / sumgt
    return(list(
      Sample.ID = sid, 
      Tool = tool,
      MD = manhattan,
      SMD = scaled.manhattan,
      sens = sens, 
      prec = prec, 
      F1 = F1,
      Combined   = 1 - scaled.manhattan + prec + sens,
      spec       = spec,
      scaled_L2 = L2 / sumgt,
      KL        = KL,
      sum_sens_spec = sens + spec
    ))
  }

  zz <- parallel::mclapply(1:nrow(xx), one.row, mc.cores = mc_cores)
  return(zz)
}

add_tool_etc <- function(exposure, tool) {
  tmp <- t(exposure)
  df1 <- data.frame(
    Sample.ID = rownames(tmp),
    Tool = tool
  )
  df2 <- cbind(df1, tmp)
  rownames(df2) <- NULL
  return(df2)
}
