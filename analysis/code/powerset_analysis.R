library(rje)
library(cosmicsig)
library(BiocParallel)
library(parallel)

source("analysis/code/reoptimize_similarities_from_exposures.R")

stopifnot(packageVersion("cosmicsig") == '1.2.0')
# To get 1.2.0:
# remotes::install_github(repo = "Rozen-Lab/cosmicsig", ref = "v1.2.0-branch")

cosmic_sig_matrix = function(version = "3.4") {
  if (version == "3.4") {
    retval = COSMIC_v3.4$signature$GRCh37$SBS96
    retval = cbind(
      retval,
      COSMIC_v3.3$signature$GRCh37$SBS96[ , "SBS40", drop = FALSE]
    )
    return(retval)
  } else if (version == "3.0") {
    return(COSMIC_v3.0$signature$GRCh37$SBS96)
  } else {
    stop("argument version has unknown value: ", version)
  }
}

if (FALSE) {
global_sbs_sig_matrix = cosmicsig::COSMIC_v3.4$signature$GRCh37$SBS96
global_sbs_sig_matrix = cbind(
  global_sbs_sig_matrix,
  cosmicsig::COSMIC_v3.3$signature$GRCh37$SBS96[ , "SBS40", drop = FALSE]
)
}

experiment = function() {
  foo = pcawg_stomach_exposures()
  foo = foo[rowSums(foo) > 0 , ]
  
  mincol = function(col) {
    mincol = min(col[col > 0])
    return(mincol / sum(col))
  }
  
  retval = apply(X = foo, MARGIN = 2, mincol)

  return(retval)
  
}

experiment2 = function() {
  mins = experiment()
  for (cut in c(0.01, 0.02, 0.03, 0.03)) {

    cat(cut, " ", sum(mins >= cut), "\n")
  }
}


explore_powerset_noprint <- function(test_spectrum,
                                     all_sig_names,
                                     cosine_threshold) {
  # The name of the test_spectrum is used as a file name,
  # and :: is not allowed
  colnames(test_spectrum) <-
    strsplit(colnames(test_spectrum), "::", fixed = TRUE)[[1]][2]
  # browser()
  message("starting ", colnames(test_spectrum))
  
  sbs_sigs_matrix <- cosmic_sig_matrix()[, all_sig_names, drop = FALSE]
  
  sig_sets <- rje::powerSet(all_sig_names)
  # Get rid of empty set:
  sig_sets <- sig_sets[-1]
  
  reconstruct_and_similarities <- function(spectrum, sig_names) {
    sig_matrix <- sbs_sigs_matrix[, sig_names, drop = FALSE]
    
    q_exp <- mSigTools::optimize_exposure_QP(spectrum, sig_matrix)
    # browser()
    tmpr = sig_matrix %*% q_exp
    q_reconstruct <- round(sig_matrix %*% q_exp)
    attributes(q_reconstruct) <- attributes(spectrum)
    colnames(q_reconstruct) <- "q_recon"
    q_cosine <- mSigAct::cossim(spectrum, q_reconstruct)
    # noround_cosine = mSigAct::cossim(spectrum, tmpr)
    return(list(
      sig_names = paste(sig_names, collapse = ","),
      sig_num = length(sig_names),
      q_exp = q_exp,
      q_reconstruct = q_reconstruct,
      q_cosine = q_cosine # ,
      # noround_cosine = noround_cosine
    ))
  }
  
  ttime = system.time({
    if (Sys.info()["sysname"] == "Windows") {
      param <- SnowParam(workers = 5, type = "SOCK")
      all_recs <- 
        bplapply(sig_sets, 
                 function(x) reconstruct_and_similarities(test_spectrum, x),
                 BPPARAM = param)
    } else {
      all_recs = parallel::mclapply(
        sig_sets, 
        function(x) reconstruct_and_similarities(test_spectrum, x),
        mc.cores = 1)
    }
  })
  # browser()
  message("finished ", colnames(test_spectrum))

  to_examine <- unlist(lapply(all_recs, `[[`, "q_cosine")) >= cosine_threshold
  if (sum(to_examine) == 0) {
    return(list(reconstructions = NA, exposures = NA))
  }
  # OLD noround_examine = unlist(lapply(all_recs, `[[`, "noround_cosine")) >= cosine_threshold
  # OLD message("round len ", sum(to_examine))
  # OLD message("noround len ", sum(noround_examine))  
  
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
  
  # mSigTools::optimize_exposure_QP sometimes attributes exposures of < 0.5
  # mutation to a signature. "no_zero_exposure" is TRUE if optimize_exposure_QP
  # did not make assignments < 0.5. Attributions with such small assignments are
  # redundant.
  to_use <- unlist(lapply(all_exp, `[[`, "no_zero_exposure"))
  q_reconstructions <- q_reconstructions[, to_use, drop = FALSE]
  all_exp <- all_exp[to_use]
  return(list(reconstructions = q_reconstructions, exposures = all_exp))
}

explore_powerset_and_print <- function(test_spectrum,
                             all_sig_names,
                             cosine_threshold,
                             xlabels,
                             dir_name,
                             plot_pdf = TRUE) {
  if (plot_pdf && !dir.exists(dir_name)) {
    dir.create(path = dir_name, recursive = FALSE)
  }
  # browser()
  tmplist <-
    explore_powerset_noprint(
      test_spectrum = test_spectrum,
      all_sig_names = all_sig_names,
      cosine_threshold = cosine_threshold
    )
  
  if (class(tmplist$exposure) == "logical") { # Then this was an NA
    return(tmplist)
  }
  # browser()
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
        "#E69F00", "#000000", "#229966"
      )
    names(color_sets) <-
      c(
        "SBS1", "SBS2", "SBS3", "SBS5", "SBS26", "SBS13", "SBS15", "SBS17a",
        "SBS17b", "SBS18", "SBS21", "SBS40", "SBS41", "SBS44", "SBS58"
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
        # border = NA,
        col = col,
        ylab = "", legend.x = 1.2
      )
    }
    
    grDevices::dev.off()
    return(num_printed)
  }

number_reconstructions <- function(yy) {
  if (class(yy) == "logical") {
    return(0) # This is an NA
  } else {
    return(length(yy))
  }
}

filter_exposure_no_low_fraction_contribution = 
  function(exposure_record, fraction) {
    if (is.logical(exposure_record)) {
      return(FALSE)
    }
    exposure = exposure_record[["q_exp_matrix"]]
    exposure_fraction = exposure / sum(exposure)
    
    stopifnot(all.equal(sum(exposure_fraction), 1.0))
    return(all(exposure_fraction > fraction))
  }


exposure_exceeds_proportion_filter = function(exposure, proportion_threshold) {
  exposure_fraction = exposure / sum(exposure)
  
  stopifnot(all.equal(sum(exposure_fraction), 1.0))
  return(all(exposure_fraction > proportion_threshold))
}
