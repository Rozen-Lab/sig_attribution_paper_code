
library(parallel)
library(mSigTools)
library(philentropy)

gather_stats_any = function(mutation_type) {
  
  syn_exp_files <- list.files(
    path = file.path("analysis/raw_output", mutation_type),
    full.names = TRUE, recursive = TRUE,
    pattern = "^inferred_exposures.csv"
  )
  tool_names <- basename(sub("/syn.*", "", syn_exp_files))
  
  output_dir_syn <- file.path("analysis/summary", mutation_type) # Where to put the summary
  
  total_cores <- parallel::detectCores()
  cores_to_use <- total_cores / 2
  
  compare_syn_results(
    dataset = mutation_type,
    syn_exp_files = syn_exp_files,
    tool_names = tool_names,
    output_dir = output_dir_syn,
    cancer_types = NULL, # Dead code?
    mc_cores = cores_to_use
  )
  
}

compare_syn_results <-
  function(dataset, syn_exp_files, tool_names, output_dir,
           data_top_folder_name = "synthetic_data",
           cancer_types = NULL, 
           mc_cores = 1) {
    
    all_exposures <- lapply(syn_exp_files, FUN = function(file) {
      return(mSigTools::read_exposure(file))
    })
    names(all_exposures) <- tool_names

    sig_universe_mask_file = 
      file.path("synthetic_data",
                dataset, # {SBS, DBS, ID}
                "sig_universe_mask.csv")
    sig_universe_mask = mSigTools::read_exposure(sig_universe_mask_file)
 
    exp_gt <- mSigTools::read_exposure(
      file = file.path(
        data_top_folder_name,
        dataset, "ground.truth.syn.exposures.csv"
      )
    )
    
    exp_gt_df <- add_tool_etc(
      exposure = exp_gt,
      tool = "Ground-Truth"
    )
    
    exp_tool_dfs <-
      lapply(names(all_exposures), 
             FUN = function(name) {
               exp_tool_df <- add_tool_etc(
                 exposure = all_exposures[[name]],
                 tool = name
               )
               return(exp_tool_df)
             })
    
    exp_df_all <- c(list(exp_gt_df), exp_tool_dfs)
    exposure_all <- do.call(dplyr::bind_rows, exp_df_all)
    exposure_all[is.na(exposure_all)] <- 0
    
    if (!is.null(cancer_types)) {
      stop("Dead code?")
      pattern <- paste(cancer_types, collapse = "|")
      indices <-
        grep(pattern = pattern, x = exposure_all$Sample.ID)
      exposure_all <- exposure_all[indices, ]
    }
    
    all_inferred_exp_path =
      file.path(output_dir, 
                paste0("all_inferred_exposures_", dataset, ".csv"))
    message("writing ", all_inferred_exp_path)
    data.table::fwrite(exposure_all, all_inferred_exp_path)

    compute_and_write_stats(
      exposure_all = exposure_all,
      sig_universe_mask = sig_universe_mask,
      output_dir = output_dir,
      mc_cores = mc_cores,
      mutation_type = dataset
    )
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

compute_and_write_stats <- 
  function(exposure_all, sig_universe_mask,
           output_dir, mc_cores, mutation_type) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  message("writing statistics from compute_and_write_stats to ", output_dir)

  stats <- all_measures(exposure_all, 
                        sig_universe_mask,
                        mc_cores = mc_cores
  ) # Use mc_cores = 1 for debugging / testing
  
  s2 <- t(matrix(unlist(stats), nrow = length(stats[[1]]), byrow = FALSE))
  colnames(s2) <- c(
    "Sample.ID", "Tool", "MD", "SMD",
    "sens", "prec", "F1", "Combined", "spec", "scaled_L2", "KL"
  )
  
  as_tibble(s2) %>%
    filter(Tool != "Ground-Truth") %>%
    mutate(
      SMD = as.numeric(SMD),
      sens = as.numeric(sens),
      prec = as.numeric(prec),
      F1 = as.numeric(F1),
      Combined = as.numeric(Combined),
      cancer.type = sub("::.*", "", Sample.ID)
    ) -> assessment_each_sample
  
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
      m.F1         = mean(F1)
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
      m.F1         = mean(F1)
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
all_measures <- 
  function(xx, # A data.frame containing expsures. 
               #Samples are in rows and signatures are in columns
          sig_universe_mask,
          mc_cores) { # Number of cores to use
  # mc_cores <- mSigAct:::Adj.mc.cores(mc_cores)
  # mc_cores = 1 # debugging

  one.row <- function(ii) {
    # ii is a row index in xx

    if (ii == 8101) {
      # browser()
      # message("row 8101")
    }
    rr <- xx[ii, ]
    # rr is a 1-row tibble
    sid <- dplyr::pull(rr, Sample.ID)
    tool <- dplyr::pull(rr, Tool)
    
    if (tool == "Ground-Truth") {
      return(list(
        Sample.ID = sid, Tool = tool,
        MD = NA, SMD = NA, sens = NA, prec = NA, F1 = NA, Combined = NA,
        spec = NA, scaled_L2 = NA, KL = NA # Added for revision of paper
      ))
    } else {
      # browser()
      gt <- dplyr::filter(  # Get the ground truth exposures for one sample
        xx, Sample.ID == sid,
        Tool == "Ground-Truth"
      )
      one_sig_universe_mask = t(sig_universe_mask[, sid, drop = FALSE] == 1)
      gt <- dplyr::mutate(gt, Sample.ID = NULL, Tool = NULL)
      me <- dplyr::mutate(rr, Sample.ID = NULL, Tool = NULL)
      stopifnot(colnames(gt) == colnames(me))
      stopifnot(colnames(gt) == colnames(one_sig_universe_mask))
      if (!all(dim(gt) == dim(me))) {
        print(dim(gt))
        print(dim(me))
        
        browser()
      }
      
      # Only look at signatures that were offered in the universe of signatures
      # for this sample.
      gt = gt[ , one_sig_universe_mask]
      me = me[ , one_sig_universe_mask]
      me = unlist(me)
      gt = unlist(gt)
      # browser()
      KL = philentropy::kullback_leibler_distance(
        P = gt / sum(gt), # The actual distribtion
        Q = me / sum(me), # Our estimate
        testNA = FALSE, unit = "log2", epsilon = 0.001)
      
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
      
      spec <- TN / N

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
        Sample.ID = sid, Tool = tool,
        MD = manhattan,
        SMD = scaled.manhattan,
        sens = sens, prec = prec, F1 = F1,
        Combined = 1 - scaled.manhattan + prec + sens,
        spec = spec,
        scaled_L2 = L2 / sumgt,
        KL = KL
      ))
    }
  }
  zz <- parallel::mclapply(1:nrow(xx), one.row, mc.cores = mc_cores)
  return(zz)
}


