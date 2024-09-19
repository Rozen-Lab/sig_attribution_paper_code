library(rje)
library(tidyverse)
library(writexl)
library(rvest)
library(dplyr)
library(stringr)
# library(ggplot2)
library(parallel)
library(gtools)
library(data.table)
library(mSigTools)
library(mSigAct)

library(ICAMS) # remotes::install_github("steverozen/ICAMS", ref = "v3.0.8-branch")

source("analysis/code/generic_analysis.R")

get_exposure <- function(exposure, tool) {
  tmp <- t(exposure)
  df1 <- data.frame(
    Sample.ID = rownames(tmp),
    Tool = tool
  )
  df2 <- cbind(df1, tmp)
  rownames(df2) <- NULL
  return(df2)
}

all_stats <- function(xx, mc_cores) {
  mc_cores <- mSigAct:::Adj.mc.cores(mc_cores)

  one.row <- function(ii) {
    # ii is a row index in xx
    
    if (ii == 8101) {
      browser()
      message("row 8101")
    }
    rr <- xx[ii, ]
    # rr is a 1-row tibble
    sid <- dplyr::pull(rr, Sample.ID)
    tool <- dplyr::pull(rr, Tool)

    if (tool == "Ground-Truth") {
      return(list(
        Sample.ID = sid, Tool = tool,
        MD = NA, SMD = NA, sens = NA, prec = NA, F1 = NA, Combined = NA
      ))
    } else {
      gt <- dplyr::filter(  # Get the ground truth exposures for one sample
        xx, Sample.ID == sid,
        Tool == "Ground-Truth"
      )
      gt <- dplyr::mutate(gt, Sample.ID = NULL, Tool = NULL)
      me <- dplyr::mutate(rr, Sample.ID = NULL, Tool = NULL)
      stopifnot(colnames(gt) == colnames(me))
      if (!all(dim(gt) == dim(me))) {
        print(dim(gt))
        print(dim(me))

        browser()
      }
      md <- sum(abs(gt - me))

      called.pos <- which(me > 0)
      real.pos <- which(gt > 0)

      TP <- length(intersect(called.pos, real.pos))

      P <- length(real.pos)

      sens <- TP / P

      called.neg <- which(me == 0)
      real.neg <- which(gt == 0)

      TN <- length(intersect(called.neg, real.neg))
      N <- length(real.neg)

      # spec <- TN / N
      # spec is not correct because the number of called.neg and the number of
      # real.neg are inflated by signatures not in the signature universe for a
      # given cancer type.

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
      scaled.manhattan <- md / sumgt
      return(list(
        Sample.ID = sid, Tool = tool,
        MD = md, SMD = scaled.manhattan,
        sens = sens, prec = prec, F1 = F1,
        Combined = 1 - scaled.manhattan + prec + sens
      ))
    }
  }
  zz <- parallel::mclapply(1:nrow(xx), one.row, mc.cores = mc_cores)
  return(zz)
}


compute_and_write_stats <- function(exposure_all, output_dir, mc_cores) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  message("writing statistics from compute_and_write_stats to ", output_dir)
  data.table::fwrite(
    exposure_all,
    file.path(output_dir, "all_inferred_exposures.csv")
  )
  stats <- all_stats(exposure_all,
    mc_cores = mc_cores
  ) # Use mc_cores = 1 for debugging / testing

  s2 <- t(matrix(unlist(stats), nrow = length(stats[[1]]), byrow = FALSE))
  colnames(s2) <- c(
    "Sample.ID", "Tool", "MD", "SMD",
    "sens", "prec", "F1", "Combined"
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

  write.csv(summary.stats2,
    file = file.path(output_dir, "all_summary_stats.csv")
  )

  data.table::fwrite(
    assessment_each_sample,
    file = file.path(output_dir, "assessment_each_sample.csv")
  )
}


compare_syn_results <-
  function(dataset, syn_exp_files, tool_names, output_dir,
           data_top_folder_name = "synthetic_data",
           cancer_types = NULL, mc_cores = 1) {
    all_exposures <- lapply(syn_exp_files, FUN = function(file) {
      return(mSigTools::read_exposure(file))
    })

    names(all_exposures) <- tool_names

    exp_gt <- mSigTools::read_exposure(
      file = file.path(
        data_top_folder_name,
        dataset, "ground.truth.syn.exposures.csv"
      )
    )

    exp_gt_df <- get_exposure(
      exposure = exp_gt,
      tool = "Ground-Truth"
    )

    exp_tool_dfs <- lapply(names(all_exposures), FUN = function(name) {
      exp_tool_df <- get_exposure(
        exposure = all_exposures[[name]],
        tool = name
      )
      return(exp_tool_df)
    })

    exp_df_all <- c(list(exp_gt_df), exp_tool_dfs)
    exposure_all <- do.call(dplyr::bind_rows, exp_df_all)
    exposure_all[is.na(exposure_all)] <- 0

    if (!is.null(cancer_types)) {
      pattern <- paste(cancer_types, collapse = "|")
      indices <-
        grep(pattern = pattern, x = exposure_all$Sample.ID)
      exposure_all <- exposure_all[indices, ]
    }

    compute_and_write_stats(
      exposure_all = exposure_all,
      output_dir = output_dir,
      mc_cores = mc_cores
    )
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

get_gt_sig_count <- function(gt_exp_file, cancer_type) {
  gt_exp <- mSigTools::read_exposure(gt_exp_file)
  samples_one_type <-
    grep(pattern = cancer_type, x = colnames(gt_exp), value = TRUE)
  gt_exp_one_type <-
    mSigAct:::RemoveZeroActivitySig(gt_exp[, samples_one_type])
  tmp <- apply(X = gt_exp_one_type, MARGIN = 1, FUN = function(counts) {
    return(length(counts[counts > 0]))
  })
  df <- data.frame(gt_sig_id = names(tmp), gt_sig_count = tmp)
  return(df)
}

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


format_number <- function(x) {
  format(round(x, 2), nsmall = 2)
}

get_med_iqr_df <- function(summary_info) {
  median <- as.numeric(format_number(summary_info["Median"]))
  iqr <- paste0(
    format_number(summary_info["1st Qu."]),
    " to ", format_number(summary_info["3rd Qu."])
  )
  mean <- as.numeric(format_number(summary_info["Mean"]))
  df <- data.frame(Median = median, IQR = iqr, Mean = mean)
  return(df)
}

create_median_iqr_table <- function(dt) {
  dt$one_minus_SMD <- 1 - dt$SMD

  tools <- get_tool_order(dt)
  metrics <-
    c("Combined", "one_minus_SMD", "prec", "sens")

  cancer_types <- unique(dt$cancer.type)

  retval <- lapply(metrics, FUN = function(metric) {
    med_iqr_info_all_type <- lapply(tools, FUN = function(tool) {
      tmp <- summary(dt[Tool == tool, get(metric)])
      df <- get_med_iqr_df(tmp)
      colnames(df) <- paste0("All_cancer_types_", colnames(df))
      rownames(df) <- tool
      return(df)
    })
    med_iqr_df_all_type <- do.call(rbind, med_iqr_info_all_type)


    med_iqr_info_by_type <- lapply(cancer_types, FUN = function(cancer_type) {
      statistics <- lapply(tools, FUN = function(tool) {
        tmp <-
          summary(dt[Tool == tool & cancer.type == cancer_type, get(metric)])
        df <- get_med_iqr_df(tmp)
        colnames(df) <- paste0(cancer_type, "_", colnames(df))
        rownames(df) <- tool
        return(df)
      })

      statistics_df <- do.call(rbind, statistics)
      return(statistics_df)
    })

    med_iqr_df_by_type <- do.call(cbind, med_iqr_info_by_type)

    metric_df <-
      data.frame(
        Metric = metric,
        Tool = tools
      )
    med_iqr_df <- cbind(metric_df, med_iqr_df_all_type, med_iqr_df_by_type)
    return(med_iqr_df)
  })

  retval2 <- do.call(rbind, retval)

  retval2[retval2 == "one_minus_SMD"] <- "1 - scaled Manhattan distance"
  retval2[retval2 == "prec"] <- "Precision"
  retval2[retval2 == "sens"] <- "Recall"

  retval3 <- data.table::as.data.table(retval2)
  retval3 <- change_tool_names(retval3)
  return(retval3)
}

explore_powerset_noprint <- function(test_spectrum,
                                     all_sig_names,
                                     cosine_threshold) {
  # The name of the test_spectrum is used as a file name,
  # and :: is not allowed
  colnames(test_spectrum) <-
    strsplit(colnames(test_spectrum), "::", fixed = TRUE)[[1]][2]

  sbs_sigs_matrix <-
    cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96[, all_sig_names, drop = FALSE]

  sig_sets <- rje::powerSet(all_sig_names)
  # Get rid of empty set:
  sig_sets <- sig_sets[-1]

  reconstruct_and_similarities <- function(spectrum, sig_names) {
    sig_matrix <- sbs_sigs_matrix[, sig_names, drop = FALSE]

    q_exp <- mSigTools::optimize_exposure_QP(spectrum, sig_matrix)

    q_reconstruct <- round(sig_matrix %*% q_exp)
    attributes(q_reconstruct) <- attributes(spectrum)
    colnames(q_reconstruct) <- "q_recon"
    q_cosine <- mSigAct::cossim(spectrum, q_reconstruct)

    return(list(
      sig_names = paste(sig_names, collapse = ","),
      sig_num = length(sig_names),
      q_exp = q_exp,
      q_reconstruct = q_reconstruct,
      q_cosine = q_cosine
    ))
  }

  all_recs <- lapply(sig_sets, function(x) reconstruct_and_similarities(test_spectrum, x))

  to_examine <- unlist(lapply(all_recs, `[[`, "q_cosine")) >= cosine_threshold
  if (sum(to_examine) == 0) {
    return(list(reconstructions = NA, exposures = NA))
  }

  subset <- all_recs[to_examine]

  new_records <- function(xx) {
    q_exp <- deparse(xx$q_exp, width.cutoff = 500)
    tmp <- paste(
      colnames(test_spectrum),
      xx$sig_names,
      xx$sig_num,
      "cosine sim =",
      format(xx$q_cosine, digits = 5)
    )

    colnames(xx$q_reconstruct) <- tmp

    return(list(
      sig_names = xx$sig_names,
      sig_num = xx$sig_num,
      q_exp = q_exp,
      q_reconstruct = xx$q_reconstruct,
      q_cosine = xx$q_cosine
    ))
  }

  subset2 <- lapply(subset, new_records)

  subset4 <-
    subset2[
      order(unlist(lapply(subset2, `[[`, "q_cosine")), decreasing = TRUE)
    ]

  extract_exposure <- function(xx) {
    ex1 <- eval(parse(text = xx$q_exp))
    keep <- ex1 > 0.5
    q_exp_matrix <- as.matrix(ex1[keep])
    no_zero_exposure <- xx$sig_num == nrow(q_exp_matrix)
    return(list(
      sig_num = xx$sig_num,
      q_exp_matrix = q_exp_matrix,
      no_zero_exposure = no_zero_exposure,
      xx$q_cosine
    ))
  }

  all_exp <- lapply(subset4, extract_exposure)

  q_reconstructions <- do.call(cbind, lapply(subset4, `[[`, "q_reconstruct"))

  # Often mSigTools::optimize_exposure_QP sometimes attributes exposures of < 0.5
  # mutation to a signature. "no_zero_exposure" is TRUE if optimize_exposure_QP
  # did not make assignments < 0.5. Attributions with such small assignments are
  # redundant.
  to_use <- unlist(lapply(all_exp, `[[`, "no_zero_exposure"))
  q_reconstructions <- q_reconstructions[, to_use, drop = FALSE]
  all_exp <- all_exp[to_use]
  return(list(reconstructions = q_reconstructions, exposures = all_exp))
}

explore_powerset <- function(test_spectrum,
                             all_sig_names,
                             cosine_threshold,
                             xlabels,
                             dir_name,
                             plot_pdf = TRUE) {
  if (!dir.exists(dir_name)) {
    dir.create(path = dir_name, recursive = FALSE)
  }

  tmplist <-
    explore_powerset_noprint(
      test_spectrum = test_spectrum,
      all_sig_names = all_sig_names,
      cosine_threshold = cosine_threshold
    )

  if (class(tmplist$exposure) == "logical") { # Then this was an NA
    return(tmplist)
  }

  # The name of the test_spectrum is used as a file name,
  # and :: is not allowed
  colnames(test_spectrum) <-
    strsplit(colnames(test_spectrum), "::", fixed = TRUE)[[1]][2]
  if (plot_pdf) {
    plot_reconstruction_pdf(
      file = file.path(
        dir_name,
        paste(colnames(test_spectrum),
          "_similar_at_",
          cosine_threshold,
          ".pdf",
          sep = ""
        )
      ),
      spectrum = test_spectrum,
      reconstructions = tmplist$reconstructions,
      exposures = tmplist$exposures,
      xlabels = xlabels
    )
  }

  return(tmplist)
}

# Print the original spectrum and each reconstruction and the associated
# signature exposures (a.k.a. signature activities)
plot_reconstruction_pdf <-
  function(file, # Name of the output file to create
           spectrum, # Spectrum that was was reconstructed
           reconstructions, # Matrix of reconstructions, one per column
           exposures, # List of exposures and other information
           xlabels = TRUE # Passed to ICAMS::PlotCatalog
  ) {
    grDevices::pdf(
      file = file,
      width = 8.2677,
      height = 11.6929, onefile = TRUE
    )

    color_sets <-
      c(
        "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
        "#44AA99", "#882255", "#999933", "#661100", "#6699CC", "#888888",
        "#E69F00", "#000000"
      )
    names(color_sets) <-
      c(
        "SBS1", "SBS2", "SBS3", "SBS5", "SBS26", "SBS13", "SBS15", "SBS17a",
        "SBS17b", "SBS18", "SBS21", "SBS40", "SBS41", "SBS44"
      )

    # Divide the plotting screen into 6 rows and 7 columns;
    # spectrum plot will use the first 6 columns,
    # exposure plot will use the 7th column
    mat <- matrix(
      c(
        rep(1, 6), 2,
        rep(3, 6), 4,
        rep(5, 6), 6,
        rep(7, 6), 8,
        rep(9, 6), 10,
        rep(11, 6), 12
      ),
      nrow = 6,
      ncol = 7,
      byrow = TRUE
    )
    nf <- layout(mat)

    # First row is for plotting the original spectrum
    par(mar = c(4, 5.5, 5.5, 1), tck = 0)
    ICAMS::PlotCatalog(spectrum, xlabels = xlabels)
    plot.new()

    num_printed <- 0
    # Plot the reconstructed spectrum and the corresponding exposure
    for (i in 1:ncol(reconstructions)) {
      num_printed <- num_printed + 1
      par(mar = c(4, 5.5, 5.5, 1))
      ICAMS::PlotCatalog(reconstructions[, i, drop = FALSE], # ,
        xlabels = xlabels
      )

      par(mar = c(2, 1, 4, 4), tck = 0)

      exps <- exposures[[i]]$q_exp_matrix

      col <- color_sets[rownames(exps)]
      mSigTools::plot_exposure(
        exposure = exps,
        plot.sample.names = FALSE,
        cex.yaxis = 0.9,
        border = NA,
        col = col,
        ylab = "", legend.x = 1.2
      )
    }

    grDevices::dev.off()
    return(num_printed)
  }

figure_2_and_sup_table_s2 <- function(sample_name, selected_sig_names, dir_name) {
  test_spectrum <- PCAWG7::spectra$PCAWG$SBS96[, sample_name, drop = FALSE]
  r_1 <- explore_powerset(test_spectrum,
    selected_sig_names,
    cosine_threshold = 0.980,
    xlabels = FALSE,
    dir_name = dir_name,
    plot_pdf = FALSE
  )
  stopifnot(length(r_1$exposures) == 123)

  to_show1 <- which(colnames(r_1$reconstructions) ==
    "SP85251 SBS1,SBS2,SBS5,SBS15,SBS18,SBS40 6 cosine sim = 0.98305")
  to_show2 <- to_show1 + 1
  to_show3 <- which(colnames(r_1$reconstructions) ==
    "SP85251 SBS1,SBS2,SBS5,SBS13,SBS15,SBS17a,SBS17b,SBS18,SBS41,SBS44 10 cosine sim = 0.98051")
  to_show4 <- to_show3 + 1

  all_to_show <- c(to_show1, to_show2, to_show3, to_show4)

  fig_reconstructions <- r_1$reconstructions[, all_to_show]
  fig_exposures <- r_1$exposures[all_to_show]

  # Plot draft figure for later manual cleanup

  plot_reconstruction_pdf(
    file = file.path(dir_name, "figure_2_panels_B-F.pdf"),
    spectrum = test_spectrum,
    reconstructions = fig_reconstructions,
    exposures = fig_exposures,
    xlabels = FALSE
  )

  # Write draft supplementary table corresponding to the figure.

  lapply(fig_exposures, `[[`, "q_exp_matrix") %>%
    lapply(t) %>%
    lapply(as_tibble) %>%
    do.call(bind_rows, .) %>%
    as.matrix() %>%
    t() %>%
    data.frame() %>%
    cbind(rownames(.), .) -> tmp_fig
  colnames(tmp_fig) <- c("signature", colnames(fig_reconstructions))
  writexl::write_xlsx(data.frame(tmp_fig),
    path = file.path(dir_name, "sup_table_s2.xlsx")
  )
}

enumeration_of_attributions <- function(spectra, selected_sig_names, dir_name) {
  retval <- lapply(
    1:ncol(spectra),
    function(x) {
      cat(x, " ")
      explore_powerset(
        test_spectrum = spectra[, x, drop = FALSE],
        all_sig_names = selected_sig_names,
        cosine_threshold = 0.980,
        xlabels = FALSE,
        dir_name = dir_name
      )
    }
  )
}

number_reconstructions <- function(yy) {
  if (class(yy) == "logical") {
    return(0)
  } else {
    return(length(yy))
  }
}

figure_2_histogram <- function(file, num_of_recon) {
  cairo_pdf(file = file)
  par(mfrow = c(2, 1), mai = c(1.02, 1.25, 0.82, 0.42))
  hist(num_of_recon,
    breaks = seq(from = 0, to = 1400, by = 50),
    xlab = "Number of attributions with cosine similarity to spectrum > 0.98",
    ylab = "Number of stomach\nadenocarcinoma spectra",
    cex.lab = 1.15,
    main = ""
  )
  dev.off()
}

get_alternative_attribution <-
  function(sample_names, spectra, selected_sig_names) {
    retval <- lapply(
      sample_names,
      function(x) {
        cat(x, "\n")
        ret <-
          list(
            sample = x,
            explore_powerset_noprint(
              test_spectrum = spectra[, x, drop = FALSE],
              all_sig_names = selected_sig_names,
              cosine_threshold = 0
            )
          )
        return(ret)
      }
    )
    return(retval)
  }

# Compare ground truth attribution to attributions in the power set
compare_to_powerset <- function(yy) {
  s96.sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96

  sampleid <- yy[[1]]
  rest <- yy[[2]]
  target <- all_syn[, sampleid, drop = FALSE]
  gtexp <- as.matrix(all_syn_exp[, ..sampleid, drop = FALSE])
  rownames(gtexp) <- unlist(all_syn_exp[, 1])
  usedsigs <- rownames(gtexp)[gtexp[, 1] > 0]
  sigsusednames <- rownames(gtexp)
  sigused <- as.matrix(s96.sigs[, sigsusednames, drop = FALSE])
  reoptimized_exp <- mSigTools::optimize_exposure_QP(target, sigused[, usedsigs])
  reoptimized_recon <- sigused[, usedsigs] %*% reoptimized_exp
  reoptimized_sim <- mSigAct::cossim(reoptimized_recon, target)

  gtnumsigs <- length(usedsigs)
  recon <- sigused %*% gtexp
  gtcossim <- mSigAct::cossim(target, recon)
  processexps <- function(zz) {
    cosine <- zz[[4]]
    signum <- zz[["sig_num"]]
    if (cosine > reoptimized_sim) {
      if (signum >= gtnumsigs) {
        return(1)
      } else {
        return(-1)
      }
    } else {
      return(0)
    }
  }
  exps <- rest[["exposures"]]
  betterthangt <- unlist(lapply(exps, processexps))
  numbettergteq <- sum(betterthangt > 0)
  numbetterless <- sum(betterthangt < 0)
  return(list(
    sample_id = sampleid,
    ground_truth_num_sigs = gtnumsigs,
    ground_truth_sigs = paste(usedsigs, collapse = " "),
    reoptimized_sim = reoptimized_sim,
    ground_truth_sim = gtcossim,
    better_sim_and_ge_sigs = numbettergteq,
    better_sim_and_fewer_sigs = numbetterless,
    total_unique_attributions_gt_thresh = length(exps)
  ))
}

get_alter_attribution_table <- function(raw_results) {
  all.morexx <- lapply(raw_results, compare_to_powerset)
  table.results <- mSigAct:::ListOfList2Tibble(all.morexx)
  table.results <-
    dplyr::mutate(table.results,
      num_all_better =
        better_sim_and_ge_sigs + better_sim_and_fewer_sigs
    )
  table.results2 <-
    dplyr::select(table.results, sample_id, reoptimized_sim, num_all_better)
  table.results3 <-
    dplyr::mutate(table.results2, sample_id = gsub("Stomach-AdenoCA::", "", sample_id,
      fixed = TRUE
    ))
  return(table.results3)
}
