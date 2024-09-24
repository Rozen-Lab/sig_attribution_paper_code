
library(parallel)
library(mSigTools)
library(philentropy)
library(ggforce)
source("analysis/code/get_all_input.R")

gather_stats_any = function(mutation_type) {
  
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
  # browser()
  # exp_df_all <- c(list(exp_gt_df), exp_tool_dfs)
  # exposure_all <- do.call(dplyr::bind_rows, exp_df_all)
  exposure_all <- do.call(dplyr::bind_rows, exp_tool_dfs)
  exposure_all[is.na(exposure_all)] <- 0
  # browser()
  output_dir <- file.path("analysis/summary", mutation_type) # Where to put the summary
  all_inferred_exp_path =
    file.path(output_dir, 
              paste0("all_inferred_exposures_", mutation_type, ".csv"))
  message("writing ", all_inferred_exp_path)
  data.table::fwrite(exposure_all, all_inferred_exp_path)
  
  compute_and_write_stats(
    exposure_all = exposure_all,
    output_dir = output_dir,
    mutation_type = mutation_type
  )
}

compute_and_write_stats <- 
  function(exposure_all,
           output_dir, 
           mutation_type) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  message("writing statistics from compute_and_write_stats to ", output_dir)

  
  stats <- all_measures(exposure_all,
                        mutation_type)

  s2 <- t(matrix(unlist(stats), nrow = length(stats[[1]]), byrow = FALSE))
  colnames(s2) <- c(
    "Sample.ID", "Tool", "MD", "SMD",
    "sens", "prec", "F1", "Combined", "spec", "scaled_L2", "KL", "multiLLH"
  )
  
  as_tibble(s2) %>%
    filter(Tool != "Ground-Truth") %>%
    mutate(
      SMD = as.numeric(SMD),
      sens = as.numeric(sens),
      prec = as.numeric(prec),
      F1 = as.numeric(F1),
      Combined = as.numeric(Combined),
      spec = as.numeric(spec),
      scaled_L2 = 1 - as.numeric(scaled_L2),
      KL        = log2(as.numeric(KL) + 1),
      multiLLH  = as.numeric(multiLLH),
      cancer.type = sub("::.*", "", Sample.ID)
    ) -> assessment_each_sample

  # Can factor out the summrize call  
  assessment_each_sample %>%
    group_by(Tool) %>%
    summarize(
      m.Combined   = mean(Combined),
      med.Combined = median(Combined),
      sd.Combined  = sd(Combined),
      m.F1         = mean(F1),
      med.F1       = median(F1),
      sd.F1        = sd(F1),
      m.SMD        = mean(SMD),
      med.SMD      = median(SMD),
      sd.SMD       = sd(SMD),
      m.sens       = mean(sens),
      med.sens     = median(sens),
      sd.sens      = sd(sens),
      m.prec       = mean(prec),
      med.prec     = median(prec),
      sd.prec      = sd(prec),
      m.F1         = mean(F1),
      m.sens       = mean(sens),
      m.scaled_L2  = mean(scaled_L2),
      m.KL          = mean(KL),
      med.KL        = median(KL),
      mean.multiLLH = mean(multiLLH),
      med.multiLLH  = median(multiLLH)
    ) ->
    summary.stats
  
  summary.stats2 <- dplyr::arrange(summary.stats, desc(m.Combined))
  
  write.csv(
    summary.stats2,
    file = file.path(output_dir, 
                     paste0("all_summary_stats_", mutation_type, ".csv")))
  
  assessment_each_sample %>%
    group_by(Tool, cancer.type) %>%
    summarize(
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
      m.F1          = mean(F1),
      m.sens        = mean(sens),
      m.scaled_L2   = mean(scaled_L2),
      m.KL          = mean(KL),
      med.KL        = median(KL),
      mean.multiLLH = mean(multiLLH),
      med.multiLLH  = median(multiLLH)
    ) ->
    summary.stats.by.cancer.type
  
  summary.stats.by.cancer.type2 <- 
    dplyr::arrange(summary.stats.by.cancer.type, desc(m.Combined))
  
  write.csv(
    summary.stats.by.cancer.type2,
    file = file.path(
      output_dir, 
      paste0("summary_stats_by_cancer_type_", mutation_type, ".csv")))

  data.table::fwrite(
    assessment_each_sample,
    file = file.path(
      output_dir, 
      paste0("assessment_each_sample_", mutation_type, ".csv")))
}


# Compute all measures for each row (sample) in table xx,
# which contains 1 exposure vector per row (sample).
all_measures <- function(xx, # A data.frame containing expsures. 
                         #Samples are in rows and signatures are in columns
                         mutation_type) { 
  
  total_cores <- parallel::detectCores()
  mc_cores <- total_cores / 2
  
  # Number of cores to use

  mc_cores = 1 # debugging
  
  # browser()
  sig_universe_mask_file = 
    file.path("synthetic_data",
              mutation_type, # {SBS, DBS, ID}
              "sig_universe_mask.csv")
  sig_universe_mask = mSigTools::read_exposure(sig_universe_mask_file)
  
  ground_truth_exposures <- get_ground_truth_exposure(mutation_type)
  all_input = get_all_input(mutation_type)
  # browser()
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
    # gt0 <- dplyr::mutate(gt0, Sample.ID = NULL, Tool = NULL)
    me0 <- dplyr::mutate(rr, Sample.ID = NULL, Tool = NULL)
    stopifnot(colnames(me0) == colnames(one_sig_universe_mask))
    if (ncol(me0) > length(gt0)) {
      # This is because two signatures, SBS24 and SBS39 did not
      # have any mutations in any tumor
      stopifnot(setequal(
        setdiff(colnames(me0), names(gt0)), c("SBS24", "SBS39")))
      gt0["SBS24"] = 0
      gt0["SBS39"] = 0
    }
    gt0 = gt0[colnames(me0)]
    stopifnot(names(gt0) == colnames(me0))
    # stopifnot(sort(colnames(gt0)) == sort(colnames(one_sig_universe_mask)))
    # xone_sig_universe_mask = one_sig_universe_mask[, names(gt0)]
    
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
    # browser()
    KL = philentropy::kullback_leibler_distance(
      P = gt / sum(gt), # The actual distribtion
      Q = me / sum(me), # Our estimate
      testNA = FALSE, unit = "log2", epsilon = 0.001)
    # browser()
    inferred_spectrum = mSigAct::ReconstructSpectrum(
      sigs          = all_input$ground_truth_sigs,
      exp           = as.matrix(me),
      use.sig.names = TRUE
    )
    multiLLH = mSigAct:::LLHSpectrumMultinom(
      spectrum = inferred_spectrum,
      expected.counts = expected_spectra[ , sid]
    )
    
    L2 <- philentropy::distance(
      rbind(gt, me), method = "euclidean", mute.message = TRUE)
    
    manhattan <- sum(abs(gt - me))
    
    called.pos <- which(me > 0)
    real.pos <- which(gt > 0)
    
    TP <- length(intersect(called.pos, real.pos))
    
    P <- length(real.pos)
    
    sens <- TP / P
    
    called.neg <- which(me == 0)
    real.neg <- which(gt == 0)
    
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
      # browser()
      # message("0 denominator")
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
      multiLLH  = multiLLH
    ))
  }
  # browser()
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


