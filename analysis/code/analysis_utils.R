library(rje)
library(tidyverse)
library(writexl)
library(rvest) # for read_htm land ...?
library(dplyr)
library(stringr)
library(parallel)
library(gtools)
library(data.table)
library(mSigTools)
library(mSigAct)

library(ICAMS) # remotes::install_github("steverozen/ICAMS", ref = "v3.0.8-branch")

source("analysis/code/generic_analysis.R")


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

DO_NOT_USE_figure_2_and_sup_table_s2 <- function(sample_name, selected_sig_names, dir_name) {
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

Do_NOT_USE_enumeration_of_attributions <- function(spectra, selected_sig_names, dir_name) {
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

if (FALSE) {
  # Compare ground truth attribution to attributions in the power set
  OLD_compare_to_powerset <- function(yy) {
    s96.sigs <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96
    
    sampleid <- yy[[1]]
    rest <- yy[[2]]
    target <- all_syn[, sampleid, drop = FALSE]
    gtexp <- as.matrix(all_syn_exp[, ..sampleid, drop = FALSE])
    rownames(gtexp) <- unlist(all_syn_exp[, 1])
    usedsigs <- rownames(gtexp)[gtexp[, 1] > 0]
    sigsusednames <- rownames(gtexp)
    browser()
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
    # browser()
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
}
