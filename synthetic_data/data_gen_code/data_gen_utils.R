library(PCAWG7)
library(SynSigGen)
library(ICAMS)
library(cosmicsig)
library(gtools)
library(ggplot2)

###########################################################################

#' Remove signatures which have zero activity from an exposure matrix
#'
#' @param exposure Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @return A matrix (or data.frame) with zero activity signatures removed.
#'
remove_zero_activity_sigs <- function(exposure) {
  return(exposure[rowSums(exposure) > 0, , drop = FALSE])
}

#' Calculate number of samples for each cancer type in an exposure matrix
#'
#' @param exposure Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @importFrom PCAWG7 SplitPCAWGMatrixByTumorType
#'
#' @return A vector of integers showing the number of samples for each cancer
#'   type from \code{exposure}.
#'
calculate_num_samples <- function(exposure) {
  exposures_cancer_types <- PCAWG7::SplitPCAWGMatrixByTumorType(exposure)
  return(sapply(exposures_cancer_types, FUN = ncol))
}

#' Generate synthetic tumors from a subtype of exposure matrix
#'
#' @param seed A random seed to use.
#'
#' @param dir The directory in which to put the output; will be created if
#'   necessary.
#'
#' @param cancer_types A vector of character strings denoting different cancer
#'   types. This function will search \code{real_exposure} for exposures from
#'   tumors matching these strings. See \code{PCAWG7::CancerTypes()} for
#'   example.
#'
#' @param samples_per_caner_type Number of synthetic tumors to create for each
#'   cancer type. If it is \strong{one} number, then generate the \strong{same}
#'   number of synthetic tumors for each \code{cancer_types}. Or if it is a
#'   \strong{vector} of numbers, then generate synthetic tumors for each
#'   \code{cancer_types} accordingly to the number specified in the vector. The
#'   length and order of \code{samples_per_caner_type} should match that in
#'   \code{cancer_types}.
#'
#' @param input_sigs A matrix of signatures.
#'
#' @param real_exposure A matrix of real exposures.
#'
#' @param distribution Probability distribution used to model exposures due to
#'   one mutational signature. Can be \code{neg.binom} which stands for negative
#'   binomial distribution.
#'
#' @param sample_prefix_name Prefix name to add to the synthetic tumors.
#'
#' @param tumor_marker_name Tumor marker name to add to the synthetic tumors.
#'   E.g. "MSI-H", "POLE".
#'
#' @param sig_params Empirical signature parameters generated using real
#'   exposures irrespective of their cancer types. If there is only \strong{one}
#'   tumor having a signature in a cancer type in \code{real_exposure}, we
#'   cannot fit the \code{distribution} to only one data point. Instead, we will
#'   use the empirical parameter \code{size} from \code{sig_params}. The
#'   mutation of the only one tumor will be used as parameter \code{mu} in
#'   negative binomial distribution.
#'
#' @importFrom SynSigGen GetSynSigParamsFromExposures GenerateListOfSigParams
#'   GenerateSyntheticTumorsFromSigParams
#'
#' @return A list of three elements that comprise the synthetic data:
#' \enumerate{
#'  \item \code{ground.truth.catalog}: Spectra catalog with rows denoting
#'  mutation types and columns denoting sample names.
#'
#'  \item \code{ground.truth.signatures}: Signatures active in
#'  \code{ground.truth.catalog}.
#'
#'  \item \code{ground.truth.exposures}: Exposures of
#'  \code{ground.truth.signatures} in \code{ground.truth.catalog}.
#' }
#'
generate_subtype_syn_tumors <-
  function(seed, dir, cancer_types, samples_per_caner_type, input_sigs,
           real_exposure, distribution, sample_prefix_name, tumor_marker_name,
           sig_params) {
    subtype_sig_params_all <-
      SynSigGen::GetSynSigParamsFromExposures(
        exposures = real_exposure,
        distribution = distribution,
        cancer.type = tumor_marker_name,
        sig.params = sig_params
      )

    subtype_sig_params <-
      SynSigGen::GenerateListOfSigParams(
        real.exposures = real_exposure,
        cancer.types = cancer_types,
        distribution = distribution,
        sig.params = subtype_sig_params_all
      )

    synthetic_tumors <-
      SynSigGen::GenerateSyntheticTumorsFromSigParams(
        seed = seed,
        dir = dir,
        cancer.types = cancer_types,
        samples.per.cancer.type = samples_per_caner_type,
        input.sigs = input_sigs,
        sig.params = subtype_sig_params,
        distribution = distribution,
        sample.prefix.name = sample_prefix_name,
        tumor.marker.name = tumor_marker_name,
      )
    return(synthetic_tumors)
  }

#' Combine several exposure matrix according to cancer types
#'
#' @param exposure1 Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @param exposure2 Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @param exposure3 Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @section Note
#' \code{exposure1} should have more cancer types compared to \code{exposure2}
#' or \code{exposure3}.
#'
#' @importFrom PCAWG7 SplitPCAWGMatrixByTumorType
#'
#' @importFrom SynSigGen MergeExposures
#'
#' @return A combined exposure matrix according to cancer type.
#'
combine_exposure <- function(exposure1, exposure2, exposure3 = NULL) {
  exposure1_cancer_types <- PCAWG7::SplitPCAWGMatrixByTumorType(exposure1)
  exposure2_cancer_types <- PCAWG7::SplitPCAWGMatrixByTumorType(exposure2)
  if (!is.null(exposure3)) {
    exposure3_cancer_types <- PCAWG7::SplitPCAWGMatrixByTumorType(exposure3)
  }

  exposure_list <-
    lapply(names(exposure1_cancer_types), FUN = function(x) {
      cancer_type <- x
      exposure1_cancer_type <- exposure1_cancer_types[[cancer_type]]
      exposure2_cancer_type <- exposure2_cancer_types[[cancer_type]]

      if (!is.null(exposure3)) {
        exposure3_cancer_type <- exposure3_cancer_types[[cancer_type]]
        list_of_exposure <- list(
          exposure1_cancer_type, exposure2_cancer_type,
          exposure3_cancer_type
        )
      } else {
        list_of_exposure <- list(exposure1_cancer_type, exposure2_cancer_type)
      }

      list_of_exposure2 <- Filter(Negate(is.null), list_of_exposure)
      merged_exposure <-
        SynSigGen::MergeExposures(list.of.exposures = list_of_exposure2)
      return(merged_exposure)
    })
  new_exposure <- SynSigGen::MergeExposures(list.of.exposures = exposure_list)
  return(new_exposure)
}

#' Write the signature parameters of real and synthetic exposures
#'
#' @inheritParams generate_subtype_syn_tumors
#'
#' @param synthetic_exposure A matrix of synthetic exposures.
#'
#' @param mutation_type Mutation type of the data, e.g. "SBS96", "ID".
#'
#' @importFrom PCAWG7 SplitPCAWGMatrixByTumorType
#'
#' @importFrom SynSigGen GenerateListOfSigParams WriteSynSigParams
#'   GetSynSigParamsFromExposures
#'
write_sig_params <-
  function(dir, real_exposure, synthetic_exposure, cancer_types, distribution,
           sig_params, sample_prefix_name, mutation_type) {
    # Write parameters into files
    froot <- file.path(dir, "parameters")
    dir.create(path = froot, showWarnings = FALSE, recursive = TRUE)

    # Getting empirical estimates of key parameters describing exposures due to
    # signatures
    params <- SynSigGen::GenerateListOfSigParams(
      real.exposures = real_exposure,
      cancer.types = cancer_types,
      distribution = distribution,
      sig.params = sig_params
    )

    syn_exposure_cancer_types <-
      PCAWG7::SplitPCAWGMatrixByTumorType(synthetic_exposure)

    lapply(cancer_types, FUN = function(x) {
      one_cancer_type <- x
      parms <- params[[one_cancer_type]]

      parm_file <- file.path(froot, paste0(
        sample_prefix_name, one_cancer_type,
        ".params.", mutation_type, ".csv"
      ))
      cat("# Original paramaters\n", file = parm_file)
      # Suppress warning on column names on append
      suppressWarnings(
        SynSigGen::WriteSynSigParams(parms, parm_file,
          append = TRUE,
          col.names = NA
        )
      )
      index <- grep(
        pattern = one_cancer_type,
        x = names(syn_exposure_cancer_types)
      )
      syn_exp <- syn_exposure_cancer_types[[index]]

      # Sanity check; we regenerate the parameters from the synthetic exposures.
      check_params <-
        SynSigGen::GetSynSigParamsFromExposures(
          exposures = syn_exp,
          distribution = distribution,
          cancer.type = one_cancer_type,
          sig.params = sig_params
        )

      # check_params should be similar to parms
      cat("# Parameters derived from synthetic exposures\n",
        file = parm_file, append = TRUE
      )

      missing_sig_names <- setdiff(colnames(parms), colnames(check_params))
      if (length(missing_sig_names) > 0) {
        check_param2 <- matrix(NA, nrow = dim(parms)[1], ncol = dim(parms)[2])
        dimnames(check_param2) <- dimnames(parms)
        check_param2[, colnames(check_params)] <- check_params
        check_params <- check_param2
      } else {
        check_params <- check_params[, colnames(parms), drop = FALSE]
      }

      suppressWarnings(
        SynSigGen::WriteSynSigParams(check_params, parm_file,
          append = TRUE, col.names = NA
        )
      )

      if (length(missing_sig_names) > 0) {
        cat("# Some signatures not represented in the synthetic data:\n",
          file = parm_file, append = TRUE
        )
        cat("#", missing_sig_names, "\n", file = parm_file, append = TRUE)
      }

      cat("# Difference between original parameters and parameters",
        "derived from synthetic exposures\n",
        file = parm_file, append = TRUE
      )
      SynSigGen::WriteSynSigParams(parms - check_params, parm_file,
        append = TRUE, col.names = NA
      )
    })
    invisible(TRUE)
  }

#' Get the individual cancer type exposure and parameters from a large exposure
#' matrix
#'
#' @inheritParams generate_subtype_syn_tumors
#'
#' @importFrom PCAWG7 SplitPCAWGMatrixByTumorType
#'
#' @importFrom SynSigGen GenerateListOfSigParams
#'
#' @return A list of two elements:
#' \enumerate{
#'  \item \code{exposures}: A list of exposures for tumors in each cancer type.
#'  \item \code{sig_params}: A list of signature parameters for tumors in each
#'  cancer type.
#' }
get_exposure_info <- function(exposure, distribution, sig_params) {
  exposure_cancer_types <- PCAWG7::SplitPCAWGMatrixByTumorType(exposure)
  cancer_types <- names(exposure_cancer_types)
  sig_params_cancer_types <-
    SynSigGen::GenerateListOfSigParams(
      real.exposures = exposure,
      cancer.types = cancer_types,
      distribution = distribution,
      sig.params = sig_params
    )
  return(list(
    exposures = exposure_cancer_types,
    sig_params = sig_params_cancer_types
  ))
}

# Sample IDs of POLE-mutated tumors in the paper by Alexandrov, Kim,
# Haradhvala, Huang et al., 'The repertoire of Mutational Signatures in Human
# Cancer'. https://doi.org/10.1038/s41586-020-1943-3
pcawg_pole_tumor_ids <-
  c(
    "ColoRect-AdenoCA::SP81312", "ColoRect-AdenoCA::SP22031",
    "ColoRect-AdenoCA::SP16886", "ColoRect-AdenoCA::SP19295",
    "ColoRect-AdenoCA::SP17905", "ColoRect-AdenoCA::SP21400",
    "ColoRect-AdenoCA::SP18946", "ColoRect-AdenoCA::SP80615",
    "Uterus-AdenoCA::SP92659"
  )

# Sample IDs of MSI-H tumors in the paper by Alexandrov, Kim, Haradhvala, Huang
# et al., 'The repertoire of Mutational Signatures in Human Cancer'.
# https://doi.org/10.1038/s41586-020-1943-3
pcawg_msi_tumor_ids <-
  c(
    "Biliary-AdenoCA::SP99325", "ColoRect-AdenoCA::SP18310",
    "ColoRect-AdenoCA::SP17172", "ColoRect-AdenoCA::SP96133",
    "ColoRect-AdenoCA::SP110242", "ColoRect-AdenoCA::SP21017",
    "ColoRect-AdenoCA::SP22383", "ColoRect-AdenoCA::SP19215",
    "ColoRect-AdenoCA::SP18121", "ColoRect-AdenoCA::SP96118",
    "Kidney-RCC::SP102897", "Liver-HCC::SP98845",
    "Liver-HCC::SP107012", "Lymph-BNHL::SP116697",
    "Ovary-AdenoCA::SP102133", "Panc-AdenoCA::SP125732",
    "Panc-AdenoCA::SP125746", "Panc-AdenoCA::SP125770",
    "Skin-Melanoma::SP124384", "Stomach-AdenoCA::SP135287",
    "Stomach-AdenoCA::SP85122", "Stomach-AdenoCA::SP84384",
    "Stomach-AdenoCA::SP84982", "Stomach-AdenoCA::SP84392",
    "Stomach-AdenoCA::SP84439", "Uterus-AdenoCA::SP95454",
    "Uterus-AdenoCA::SP94741", "Uterus-AdenoCA::SP90209",
    "Uterus-AdenoCA::SP90989", "Uterus-AdenoCA::SP93540",
    "Uterus-AdenoCA::SP93227", "Uterus-AdenoCA::SP94348",
    "Uterus-AdenoCA::SP94933", "Uterus-AdenoCA::SP92364",
    "Uterus-AdenoCA::SP94917", "Uterus-AdenoCA::SP92460",
    "Uterus-AdenoCA::SP89909", "Kidney-ChRCC::SP123975"
  )

#' Plot the distributions of several exposures
#'
#' @param real_exposure A matrix of real exposures.
#'
#' @param synthetic_exposure A matrix of synthetic exposures.
#'
#' @param noisy_exposure A matrix of synthetic exposures with less noise
#'   added.
#'
#' @param size Negative binomial size parameter used for adding noise to
#'   \code{noisy_exposure}.
#'
#' @param distribution Probability distribution used to model exposures due to
#'   one mutational signature. Can be \code{neg.binom} which stands for negative
#'   binomial distribution.
#'
#' @param sig_params Empirical signature parameters generated using real
#'   exposures irrespective of their cancer types. If there is only \strong{one}
#'   tumor having a signature in a cancer type in \code{real_exposure}, we
#'   cannot fit the \code{distribution} to only one data point. Instead, we will
#'   use the empirical parameter \code{size} from \code{sig_params}. The
#'   mutation of the only one tumor will be used as parameter \code{mu} in
#'   negative binomial distribution.
#'
#' @param sample_prefix_name Prefix name to add to the synthetic tumors.
plot_exposure_distribution <-
  function(real_exposure, synthetic_exposure,
           noisy_exposure, size,
           distribution,
           sig_params, sample_prefix_name) {
    real_exposure_info <-
      get_exposure_info(
        exposure = real_exposure,
        distribution = distribution,
        sig_params = sig_params
      )
    synthetic_exposure_info <-
      get_exposure_info(
        exposure = synthetic_exposure,
        distribution = distribution,
        sig_params = sig_params
      )
    noisy_exposure_info <-
      get_exposure_info(
        exposure = noisy_exposure,
        distribution = distribution,
        sig_params = sig_params
      )

    cancer_types <- names(real_exposure_info$exposures)

    retval <- sapply(cancer_types, FUN = function(x) {
      one_caner_type <- x
      one_caner_type_with_prefix <- paste0(sample_prefix_name, one_caner_type)
      real_exposure_one_type <- real_exposure_info$exposures[[one_caner_type]]
      synthetic_exposure_one_type <-
        synthetic_exposure_info$exposures[[one_caner_type_with_prefix]]
      noisy_exposure_one_type <-
        noisy_exposure_info$exposures[[one_caner_type_with_prefix]]

      params_real_one_type <- real_exposure_info$sig_params[[one_caner_type]]
      params_synthetic_one_type <-
        synthetic_exposure_info$sig_params[[one_caner_type_with_prefix]]
      params_noisy_one_type <-
        noisy_exposure_info$sig_params[[one_caner_type_with_prefix]]

      plot_exposure_dist_one_type(
        real_exposure = real_exposure_one_type,
        synthetic_exposure = synthetic_exposure_one_type,
        params_real = params_real_one_type,
        params_synthetic = params_synthetic_one_type,
        noisy_exposure = noisy_exposure_one_type,
        params_noisy = params_noisy_one_type,
        size = size,
        cancer_type = one_caner_type
      )
    })
  }

#' Get the size parameter information
#'
#' @param size The size parameter used in negative binomial distribution for
#'   adding noise.
#'
#' @return A character string showing the size parameter information.
get_size_message <- function(size) {
  if (is.null(size)) {
    return("")
  } else {
    return(paste0("\n neg.binom.size = ", size))
  }
}

#' Plot the histogram of mutation counts attributed to one mutational signature
#'
#' @param counts A vector of mutation counts attributed to one signature across
#'   length(counts) samples.
#'
#' @param title The first line of title on top of the histogram.
#'
#' @param cex_main The magnification to be used for main titles relative to the
#'   current setting of cex.
#'
#' @param params A matrix of parameters of exposures attributed to mutational
#'   signatures. The rows are
#'   \describe{
#'      \item{prob}{The proportion of tumors with the signature.}
#'
#'      \item{size}{Dispersion parameter of negative binomial distribution
#'      fitted to the exposures.}
#'
#'      \item{mu}{Mean of negative binomial distribution fitted to the
#'      exposures.}
#'   }
#'   See \code{SynSigGen::signature.params$SBS96} for example.
#'
#' @param sig_id Signature ID, e.g. "SBS1".
#'
#' @param size The size parameter used in negative binomial distribution for
#'   adding noise.
#'
draw_histogram <-
  function(counts, title, cex_main, params, sig_id, size = NULL) {
    size_message <- get_size_message(size)

    hist(counts,
      breaks = 1000,
      main = paste0(
        title,
        "\n N = ", length(counts),
        " prob = ", round(params[1, sig_id], 4), " ",
        size_message,
        "\n mu = ", round(params[3, sig_id], 2),
        "\n size = ", round(params[2, sig_id], 2)
      ),
      cex.main = cex_main,
      xlab = "counts",
      probability = TRUE
    )

    if (length(counts) > 1) {
      lines(density(counts), col = "red")
    }
  }

#' Draw kernel density plot of mutation counts attributed to one mutational
#' signature
#'
#' @inheritParams draw_histogram
#'
#' @param num_sample Number of samples to be shown at the title.
#'
#' @param xlab Label for the x axis.
#'
draw_density <- function(counts, title, cex_main, params, sig_id,
                         num_sample = length(counts),
                         size = NULL, xlab = NULL) {
  size_message <- get_size_message(size)

  plot(density(counts),
    col = "red",
    main = paste0(
      title,
      "\n N = ", num_sample,
      " prob = ", round(params[1, sig_id], 4), " ",
      size_message,
      "\n mu = ", round(params[3, sig_id], 2),
      "\n size = ", round(params[2, sig_id], 2)
    ),
    cex.main = cex_main,
    xlab = xlab
  )
}

#' Draw two kernel density plots of mutation counts attributed to mutational
#' signature
#'
#' @inheritParams draw_density
#'
#' @param counts1 A vector of mutation counts attributed to one signature across
#'   \code{length(counts1)} samples.
#'
#' @param counts2 A vector of mutation counts attributed to one signature across
#'   \code{length(counts2)} samples.
#'
#' @param legend A character or expression vector of length = 1 to appear in the
#'   legend.
#'
draw_two_density <-
  function(counts1, counts2, title, cex_main, params, sig_id,
           legend, size = NULL,
           xlab = "Two-sample Kolmogorov-Smirnov test") {
    draw_density(
      counts = counts1, title = title, cex_main = cex_main,
      params = params, sig_id = sig_id, num_sample = length(counts2),
      size = size, xlab = xlab
    )

    if (length(counts2) > 1) {
      lines(density(counts2),
        col = "blue",
        xlim = c(0, max(counts1))
      )
    }

    retval <-
      suppressWarnings(ks.test(
        x = counts1,
        y = counts2
      ))
    p_value <- round(retval$p.value, 3)
    legend("topright",
      title = paste0("p-value = ", p_value),
      legend = legend,
      col = c("red", "blue"),
      fill = c("red", "blue"),
      border = "white", bty = "n"
    )
  }

#' Plot the distributions of exposures attributed to mutational signatures in
#' one cancer type
#'
#' @param real_exposure A matrix of real exposures.
#'
#' @param synthetic_exposure A matrix of synthetic exposures.
#'
#' @param params_real Signature parameters of \code{real_exposure}.
#'
#' @param params_synthetic Signature parameters of \code{synthetic_exposure}.
#'
#' @param noisy_exposure A matrix of synthetic exposures with noise added.
#'
#' @param params_noisy Signature parameters of \code{noisy_exposure}.
#'
#' @param size1 Negative binomial size parameter used for adding noise to
#'   \code{noisy_exposure}.
#'
#' @param cancer_type A character string denoting one cancer type. See
#'   \code{PCAWG7::CancerTypes()} for examples.
#'
plot_exposure_dist_one_type <-
  function(real_exposure, synthetic_exposure,
           params_real, params_synthetic,
           noisy_exposure, params_noisy, size1,
           cancer_type) {
    real_exposure <- remove_zero_activity_sigs(real_exposure)
    synthetic_exposure <- remove_zero_activity_sigs(synthetic_exposure)
    noisy_exposure <- remove_zero_activity_sigs(noisy_exposure)

    for (i in rownames(synthetic_exposure)) {
      real_exposure_non_zero <- real_exposure[i, which(real_exposure[i, ] > 0)]
      synthetic_exposure_non_zero <-
        synthetic_exposure[i, which(synthetic_exposure[i, ] > 0)]
      if (i %in% rownames(noisy_exposure)) {
        noisy_exposure_non_zero <-
          noisy_exposure[i, which(noisy_exposure[i, ] > 0)]
      } else {
        noisy_exposure_non_zero <- NULL
      }

      draw_histogram(
        counts = real_exposure_non_zero,
        title = paste0(cancer_type, ".", i, ".real.exposure"),
        cex_main = 0.9,
        params = params_real,
        sig_id = i
      )

      draw_histogram(
        counts = synthetic_exposure_non_zero,
        title = paste0(cancer_type, ".", i, ".synthetic.exposure"),
        cex_main = 0.9,
        params = params_synthetic,
        sig_id = i
      )

      if (!is.null(noisy_exposure_non_zero)) {
        draw_histogram(
          counts = noisy_exposure_non_zero,
          title = paste0(cancer_type, ".", i, ".noisy.exposure"),
          cex_main = 0.85,
          params = params_noisy,
          sig_id = i,
          size = size1
        )
      } else {
        plot.new()
      }

      if (length(real_exposure_non_zero) > 1) {
        draw_density(
          counts = real_exposure_non_zero,
          title = paste0(cancer_type, ".", i, ".real.exposure"),
          cex_main = 0.9, params = params_real, sig_id = i,
        )
      } else {
        plot.new()
      }

      # Plot several densities in one graph
      if (length(real_exposure_non_zero) > 1) {
        draw_two_density(
          counts1 = real_exposure_non_zero,
          counts2 = synthetic_exposure_non_zero,
          title = paste0(
            cancer_type, ".", i, ".synthetic.exposure"
          ),
          cex_main = 0.9,
          params = params_synthetic,
          sig_id = i,
          legend = c("real.exposure", "synthetic.exposure")
        )
      } else {
        plot.new()
      }

      if (length(real_exposure_non_zero) > 1) {
        draw_two_density(
          counts1 = real_exposure_non_zero,
          counts2 = noisy_exposure_non_zero,
          title = paste0(
            cancer_type, ".", i, ".noisy.exposure"
          ),
          cex_main = 0.85,
          params = params_noisy,
          sig_id = i,
          legend = c("real.exposure", "noisy.exposure"),
          size = size1
        )
      } else {
        plot.new()
      }
    }
  }

#' Plot catalog to a PDF file
#'
#' @param catalog A catalog as defined in \code{\link{ICAMS}} with attributes
#'   added. Use \code{ICAMS::ReadCatalog} to read in a catalog.
#'
#' @param identifier Character string showing the identity of \code{catalog},
#'   e.g. "no noise", "more noise", "less noise".
#'
#' @param file The name of the PDF file to be produced.
#'
#' @param grid A logical value indicating whether to draw grid lines.
#'
#' @param upper A logical value indicating whether to draw horizontal lines and
#'   the names of major mutation class on top of graph.
#'
#' @param xlabels A logical value indicating whether to draw x axis labels.
#'
plot_catalog_to_pdf <- function(catalog, identifier, file,
                                grid = FALSE, upper = FALSE, xlabels = FALSE) {
  colnames(catalog) <- paste0(
    colnames(catalog), " (", identifier, ") ", "(counts = ",
    colSums(catalog), ")"
  )

  if (inherits(catalog, "SBS96Catalog")) {
    # Suppress plotting  x axis tick marks for SBS96Catalog
    old_par <- par(tck = 0)
    on.exit(par(old_par))
  }

  ICAMS::PlotCatalogToPdf(
    catalog = catalog,
    file = file,
    grid = grid,
    upper = upper,
    xlabels = xlabels
  )
}

#' Get cancer type information from sample IDs of tumors
#'
#' @param sample_ids Character vector of sample IDs from tumors (e.g.
#'   "Breast-AdenoCA::SP117975", "SP.Syn.Breast-AdenoCA::S.1").
#'
#' @param syn_tumor Logical value indicating whether \code{sample_ids} are from
#'   synthetic tumors.
#'
#' @return Character vector of cancer types of tumors.
#'
get_cancer_type <- function(sample_ids, syn_tumor = FALSE) {
  cancer_types <- sapply(sample_ids, FUN = function(sample_id) {
    tmp <- strsplit(sample_id, "::")[[1]][1]
    if (syn_tumor) {
      # Now the sample names for synthetic tumors do not have prefix SP.Syn.
      return(tmp)
    } else {
      return(tmp)
    }
  })
  return(cancer_types)
}

#' Calculate distances between spectra and reconstructed
#' spectra of tumors
#'
#' @param spectra Spectra catalog as a numerical matrix with rows denoting
#'   mutation types and columns denoting sample names.
#'
#' @param exposure Exposures as a numerical matrix (or data.frame) with
#'   signatures in rows and samples in columns. Rownames are taken as the
#'   signature names and column names are taken as the sample IDs.
#'
#' @param sigs A matrix of signatures.
#'
#' @param group Character string showing the group information for
#'   \code{spectra}.
#'
#' @param syn_tumor Logical value indicating whether \code{spectra} are from
#'   synthetic tumors.
#'
#' @return A data frame with row names being the sample names in \code{spectra}.
#'   The columns show reconstruction distance information between \code{spectra} and
#'   reconstructed spectra using \code{sigs} and \code{exposure}. The column
#'   names are:
#'   * euclidean: Euclidean distance
#'   * mutations: Mutations of \code{spectra}
#'   * scaled.euclidean: euclidean / mutations
#'   * cosine: Cosine similarity
#'   * group: Group information for \code{spectra}
#'   * cancer.type: Cancer type for each sample in \code{spectra}
#'
get_distance <- function(spectra, exposure, sigs, group, syn_tumor = FALSE) {
  distance <- mSigAct:::CalculateDistance(
    spectra = spectra,
    exposure = exposure,
    sigs = sigs
  )
  distance <-
    distance[, c("euclidean", "mutations", "scaled.euclidean", "cosine")]
  distance$group <- group
  distance$cancer.type <-
    get_cancer_type(sample_ids = rownames(distance), syn_tumor = syn_tumor)
  return(distance)
}

#' Read in synthetic data from a directory
#'
#' @param dir Directory path to read in the synthetic data.
#'
#' @return A list of two elements that comprise the synthetic data:
#' \enumerate{
#'  \item \code{spectra}: Spectra catalog with rows denoting
#'  mutation types and columns denoting sample names.
#'
#'  \item \code{exposure}: Exposures of samples in \code{spectra}.
#'
#' }
#'
get_syn_data_info <- function(dir) {
  syn_spectra <-
    ICAMS::ReadCatalog(file = file.path(
      dir,
      "ground.truth.syn.catalog.csv"
    ))
  syn_exposure <-
    mSigTools::read_exposure(file = file.path(
      dir,
      "ground.truth.syn.exposures.csv"
    ))
  return(list(
    spectra = syn_spectra, exposure = syn_exposure
  ))
}

#' Add noise to noiseless synthetic data with different negative-binomial size
#' parameter
#'
#' @param seed A random seed to use.
#'
#' @param exposure Noiseless synthetic exposure.
#'
#' @param sigs A matrix of signatures.
#'
#' @param nb_sizes Character vector of numerical values specifying different
#'   negative-binomial size parameter.
#'
#' @param mc_cores The number of cores to use.
#'
#' @return A list with the same length as \code{nb_sizes}. Each element in the
#'   list is a list which contains three elements that comprise the synthetic
#'   data with noise added:
#'   \enumerate{
#'    \item \code{spectra}: Spectra catalog with noise added.
#'
#'    \item \code{exposure}: Exposures of \code{sigs} in \code{spectra}.
#'
#'    \item \code{sigs}: Signatures active in \code{spectra}.
#' }
#'
generate_noisy_data <- function(seed, exposure, sigs, nb_sizes, mc_cores = 1) {
  retval <- parallel::mclapply(nb_sizes, FUN = function(size) {
    output_dir <- file.path(tempdir(), size)
    tmp <-
      SynSigGen::GenerateNoisyTumors(
        seed = seed,
        dir = output_dir,
        input.exposure = exposure,
        signatures = sigs,
        n.binom.size = size
      )
    # Delete temp dir
    unlink(x = output_dir, recursive = TRUE)
    noise_data <- list(
      spectra = tmp$spectra, exposure = tmp$exposures,
      sigs = sigs
    )
    return(noise_data)
  }, mc.cores = mc_cores)
  names(retval) <- paste0("size_", nb_sizes)
  return(retval)
}

#' Calculate reconstruction distances for a list of synthetic data
#'
#' @param list_of_syn_data A list of synthetic data.
#'
#' @param mc_cores The number of cores to use.
#'
#' @return A list of data frames showing the scaled Euclidean distance
#'   information. See the return value of function \code{get_distance}
#'   for more details.
#'
get_multiple_syn_distances <- function(list_of_syn_data, mc_cores = 1) {
  retval <-
    parallel::mclapply(seq_len(length(list_of_syn_data)), FUN = function(index) {
      one_syn_data <- list_of_syn_data[[index]]
      group <- names(list_of_syn_data)[index]
      out <- get_distance(
        spectra = one_syn_data$spectra,
        exposure = one_syn_data$exposure,
        sigs = one_syn_data$sigs,
        group = group,
        syn_tumor = TRUE
      )
      return(out)
    }, mc.cores = mc_cores)
  names(retval) <- names(list_of_syn_data)
  return(retval)
}

#' Create one violin plot comparing different groups
#'
#' @param distance_df A data frame which has the reconstruction distance
#'   information from multiple groups.
#'
#' @param title Title of the violin plot.
#'
#' @param ylim Y-axis limit of the violin plot.
#'
#' @return A list of one ggplot object.
#'
one_violinplot <- function(distance_df, title, ylim, p_label = "p.signif") {
  # Customize x-axis tick labels
  x_lables <- unique(distance_df$group)
  x_lables <-
    gsub(pattern = "real", replacement = "Actual\ndata", x = x_lables)
  x_lables <-
    gsub(pattern = "size_", replacement = "Size=", x = x_lables)

  safe_colorblind_palette <-
    c(
      "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
      "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"
    )

  plot_object <-
    ggpubr::ggviolin(
      data = distance_df,
      x = "group", y = "cosine", color = "group", fill = "group",
      palette = safe_colorblind_palette
    ) + geom_boxplot(width = 0.1, color = "grey") +
    ggpubr::stat_compare_means(
      label = p_label, ref.group = "real", method = "wilcox.test",
      vjust = -2.2
    ) +
    ggplot2::ggtitle(title) + ggplot2::ylab("Cosine similarity") +
    ggplot2::theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    ggplot2::scale_x_discrete(labels = x_lables)

  return(list(plot_object))
}

#' Create multiple violin plots comparing different groups
#'
#' @param distance_df A data frame which has the reconstruction distance
#'   information from multiple groups.
#'
#' @param data_type Mutation type of the data which generated \code{distance_df}
#'   (e.g. "SBS", "indel").
#'
#' @param ylim Y-axis limit of the violin plot.
#'
#' @return A list of multiple ggplot objects.
#'
create_violinplots <- 
  function(distance_df, data_type, ylim, p_label = "p.signif", 
           plot_cancer_type = FALSE) {
  cancer_types <- unique(distance_df$cancer.type)
  # First create violin plot for all cancer types
  plot_objects <-
    one_violinplot(
      distance_df = distance_df,
      title = paste0("All cancer types (", data_type, ")"),
      ylim = ylim,
      p_label = p_label
    )
  
  if (plot_cancer_type) {
    # Create violin plots for individual cancer type
    for (cancer_type in cancer_types) {
      distance_df_one_type <-
        distance_df[distance_df$cancer.type == cancer_type, ]
      tmp <-
        one_violinplot(
          distance_df = distance_df_one_type,
          title = paste0(cancer_type, " (", data_type, ")"),
          ylim = ylim,
          p_label = p_label
        )
      plot_objects <- c(plot_objects, tmp)
    }
  }

  return(plot_objects)
}

#' Plot multiple ggplot objects to one PDF
#'
#' @param plot_objects A list of ggplot objects.
#'
#' @param file File path of the PDF to be created.
#'
#' @param nrow Number of rows of plot objects in one page.
#'
#' @param ncol Number of columns of plot objects in one page.
#'
#' @param width Width of the plot size.
#'
#' @param height Height of the plot size.
#'
#' @param units Units("in", "cm", "mm", or "px") of the plot size.
#'
ggplot_to_pdf <-
  function(plot_objects, file, nrow, ncol, width, height, units) {
    ggplot2::ggsave(
      filename = file,
      plot = gridExtra::marrangeGrob(plot_objects, nrow = nrow, ncol = ncol),
      width = width, height = height, units = units
    )
  }

tcga_pole_tumor_ids <- c(
  "Breast-cancer::TCGA-AN-A046-01A-21W-A050-09",
  "ColoRect-AdenoCa::TCGA-AA-3510-01A-01D-1408-10",
  "ColoRect-AdenoCa::TCGA-AA-3984-01A-02D-1981-10",
  "ColoRect-AdenoCa::TCGA-AA-A010-01A-01D-A17O-10",
  "ColoRect-AdenoCa::TCGA-AZ-4315-01A-01D-1408-10",
  "ColoRect-AdenoCa::TCGA-AG-3892-01A-01W-1073-09",
  "ColoRect-AdenoCa::TCGA-AG-A002-01A-01W-A00K-09",
  "ColoRect-AdenoCa::TCGA-EI-6917-01A-11D-1924-10",
  "Stomach-AdenoCa::TCGA-BR-8680-01A-11D-2394-08"
)

tcga_msi_tumor_ids <- c(
  "Breast-cancer::TCGA-A8-A09Z-01A-11W-A019-09",
  "Breast-cancer::TCGA-AN-A0AK-01A-21W-A019-09",
  "Breast-cancer::TCGA-AO-A128-01A-11D-A10M-09",
  "Breast-cancer::TCGA-BH-A0HA-01A-11D-A12Q-09",
  "Breast-cancer::TCGA-BH-A18G-01A-11D-A12B-09",
  "Breast-cancer::TCGA-D8-A1J8-01A-11D-A13L-09",
  "Breast-cancer::TCGA-D8-A1XK-01A-21D-A14K-09",
  "Breast-cancer::TCGA-D8-A1XQ-01A-11D-A14K-09",
  "Breast-cancer::TCGA-EW-A1IZ-01A-11D-A188-09",
  "ColoRect-AdenoCa::TCGA-5M-AAT6-01A-11D-A40P-10",
  "ColoRect-AdenoCa::TCGA-A6-2686-01A-01D-1408-10",
  "ColoRect-AdenoCa::TCGA-A6-3809-01A-01D-A270-10",
  "ColoRect-AdenoCa::TCGA-A6-5661-01A-01D-1650-10",
  "ColoRect-AdenoCa::TCGA-A6-5665-01A-01D-1650-10",
  "ColoRect-AdenoCa::TCGA-A6-6653-01A-11D-1771-10",
  "ColoRect-AdenoCa::TCGA-A6-6780-01A-11D-A270-10",
  "ColoRect-AdenoCa::TCGA-A6-6781-01A-22D-A270-10",
  "ColoRect-AdenoCa::TCGA-AA-3492-01A-01D-1408-10",
  "ColoRect-AdenoCa::TCGA-AA-3663-01A-01D-1719-10",
  "ColoRect-AdenoCa::TCGA-AA-3713-01A-21D-1719-10",
  "ColoRect-AdenoCa::TCGA-AA-3715-01A-01W-0900-09",
  "ColoRect-AdenoCa::TCGA-AA-3811-01A-01W-0995-10",
  "ColoRect-AdenoCa::TCGA-AA-3815-01A-01W-0995-10",
  "ColoRect-AdenoCa::TCGA-AA-3821-01A-01W-0995-10",
  "ColoRect-AdenoCa::TCGA-AA-3833-01A-01W-0900-09",
  "ColoRect-AdenoCa::TCGA-AA-3845-01A-01W-0995-10",
  "ColoRect-AdenoCa::TCGA-AA-3864-01A-01D-1981-10",
  "ColoRect-AdenoCa::TCGA-AA-3877-01A-01W-0995-10",
  "ColoRect-AdenoCa::TCGA-AA-3947-01A-01D-1981-10",
  "ColoRect-AdenoCa::TCGA-AA-3949-01A-01W-0995-10",
  "ColoRect-AdenoCa::TCGA-AA-3950-01A-02W-0995-10",
  "ColoRect-AdenoCa::TCGA-AA-3966-01A-01W-1073-09",
  "ColoRect-AdenoCa::TCGA-AA-A01P-01A-21D-A17O-10",
  "ColoRect-AdenoCa::TCGA-AA-A01R-01A-21D-A17O-10",
  "ColoRect-AdenoCa::TCGA-AA-A022-01A-21D-A17O-10",
  "ColoRect-AdenoCa::TCGA-AA-A02R-01A-01W-A00E-09",
  "ColoRect-AdenoCa::TCGA-AD-5900-01A-11D-1650-10",
  "ColoRect-AdenoCa::TCGA-AD-6889-01A-11D-1924-10",
  "ColoRect-AdenoCa::TCGA-AD-6895-01A-11D-1924-10",
  "ColoRect-AdenoCa::TCGA-AD-A5EJ-01A-11D-A28G-10",
  "ColoRect-AdenoCa::TCGA-AM-5821-01A-01D-1650-10",
  "ColoRect-AdenoCa::TCGA-AU-6004-01A-11D-1719-10",
  "ColoRect-AdenoCa::TCGA-AY-6197-01A-11D-1719-10",
  "ColoRect-AdenoCa::TCGA-AZ-4615-01A-01D-1408-10",
  "ColoRect-AdenoCa::TCGA-AZ-6598-01A-11D-1771-10",
  "ColoRect-AdenoCa::TCGA-CK-4951-01A-01D-1408-10",
  "ColoRect-AdenoCa::TCGA-CK-5913-01A-11D-1650-10",
  "ColoRect-AdenoCa::TCGA-CK-5916-01A-11D-1650-10",
  "ColoRect-AdenoCa::TCGA-CK-6746-01A-11D-1835-10",
  "ColoRect-AdenoCa::TCGA-CM-4743-01A-01D-1719-10",
  "ColoRect-AdenoCa::TCGA-CM-4746-01A-01D-1408-10",
  "ColoRect-AdenoCa::TCGA-CM-5861-01A-01D-1650-10",
  "ColoRect-AdenoCa::TCGA-CM-6162-01A-11D-1650-10",
  "ColoRect-AdenoCa::TCGA-CM-6171-01A-11D-1650-10",
  "ColoRect-AdenoCa::TCGA-CM-6674-01A-11D-1835-10",
  "ColoRect-AdenoCa::TCGA-D5-6530-01A-11D-1719-10",
  "ColoRect-AdenoCa::TCGA-D5-6540-01A-11D-1719-10",
  "ColoRect-AdenoCa::TCGA-D5-6927-01A-21D-1924-10",
  "ColoRect-AdenoCa::TCGA-D5-6928-01A-11D-1924-10",
  "ColoRect-AdenoCa::TCGA-D5-6930-01A-11D-1924-10",
  "ColoRect-AdenoCa::TCGA-DM-A1HB-01A-21D-A183-10",
  "ColoRect-AdenoCa::TCGA-F4-6570-01A-11D-1771-10",
  "ColoRect-AdenoCa::TCGA-F4-6703-01A-11D-1835-10",
  "ColoRect-AdenoCa::TCGA-F4-6856-01A-11D-1924-10",
  "ColoRect-AdenoCa::TCGA-G4-6302-01A-11D-1719-10",
  "ColoRect-AdenoCa::TCGA-G4-6304-01A-11D-1924-10",
  "ColoRect-AdenoCa::TCGA-G4-6309-01A-21D-1835-10",
  "ColoRect-AdenoCa::TCGA-G4-6320-01A-11D-1719-10",
  "ColoRect-AdenoCa::TCGA-G4-6586-01A-11D-1771-10",
  "ColoRect-AdenoCa::TCGA-G4-6588-01A-11D-1771-10",
  "ColoRect-AdenoCa::TCGA-G4-6628-01A-11D-1835-10",
  "ColoRect-AdenoCa::TCGA-NH-A5IV-01A-42D-A36X-10",
  "ColoRect-AdenoCa::TCGA-QG-A5Z2-01A-11D-A28G-10",
  "ColoRect-AdenoCa::TCGA-WS-AB45-01A-11D-A40P-10",
  "ColoRect-AdenoCa::TCGA-AG-3891-01A-01W-1073-09",
  "ColoRect-AdenoCa::TCGA-AG-A02N-01A-11W-A096-10",
  "ColoRect-AdenoCa::TCGA-EI-6507-01A-11D-1733-10",
  "ColoRect-AdenoCa::TCGA-EI-6882-01A-11D-1924-10",
  "Eso-AdenoCa::TCGA-L5-A43J-01A-12D-A247-09",
  "Eso-AdenoCa::TCGA-L5-A4OI-01A-11D-A27G-09",
  "Eso-AdenoCa::TCGA-L5-A8NM-01A-11D-A37C-09",
  "Kidney-RCC::TCGA-B0-5098-01A-01D-1421-08",
  "Liver-HCC::TCGA-4R-AA8I-01A-11D-A382-10",
  "Liver-HCC::TCGA-WQ-A9G7-01A-11D-A36X-10",
  "Ovary-AdenoCa::TCGA-59-2349-01A-01W-0799-08",
  "Skin-Melanoma::TCGA-D3-A8GM-06A-11D-A372-08",
  "Stomach-AdenoCa::TCGA-B7-5816-01A-21D-1600-08",
  "Stomach-AdenoCa::TCGA-B7-A5TI-01A-11D-A31L-08",
  "Stomach-AdenoCa::TCGA-BR-4184-01A-01D-1126-08",
  "Stomach-AdenoCa::TCGA-BR-4201-01A-01D-1126-08",
  "Stomach-AdenoCa::TCGA-BR-4256-01A-01D-1126-08",
  "Stomach-AdenoCa::TCGA-BR-4257-01A-01D-1126-08",
  "Stomach-AdenoCa::TCGA-BR-4292-01A-01D-1126-08",
  "Stomach-AdenoCa::TCGA-BR-4362-01A-01D-1158-08",
  "Stomach-AdenoCa::TCGA-BR-4363-01A-01D-1158-08",
  "Stomach-AdenoCa::TCGA-BR-4368-01A-01D-1158-08",
  "Stomach-AdenoCa::TCGA-BR-4370-01A-01D-1158-08",
  "Stomach-AdenoCa::TCGA-BR-6566-01A-11D-1800-08",
  "Stomach-AdenoCa::TCGA-BR-6802-01A-11D-1882-08",
  "Stomach-AdenoCa::TCGA-BR-6852-01A-11D-1882-08",
  "Stomach-AdenoCa::TCGA-BR-7703-01A-11D-2053-08",
  "Stomach-AdenoCa::TCGA-BR-7707-01A-11D-2053-08",
  "Stomach-AdenoCa::TCGA-BR-7851-01A-11D-2201-08",
  "Stomach-AdenoCa::TCGA-BR-8059-01A-11D-2340-08",
  "Stomach-AdenoCa::TCGA-BR-8078-01A-11D-2340-08",
  "Stomach-AdenoCa::TCGA-BR-8081-01A-11D-2340-08",
  "Stomach-AdenoCa::TCGA-BR-8284-01A-11D-2340-08",
  "Stomach-AdenoCa::TCGA-BR-8360-01A-11D-2340-08",
  "Stomach-AdenoCa::TCGA-BR-8361-01A-11D-2340-08",
  "Stomach-AdenoCa::TCGA-BR-8363-01A-11D-2340-08",
  "Stomach-AdenoCa::TCGA-BR-8368-01A-11D-2340-08",
  "Stomach-AdenoCa::TCGA-BR-8382-01A-11D-2394-08",
  "Stomach-AdenoCa::TCGA-BR-8487-01A-11D-2394-08",
  "Stomach-AdenoCa::TCGA-BR-A4QL-01A-31D-A25D-08",
  "Stomach-AdenoCa::TCGA-CD-8536-01A-11D-2340-08",
  "Stomach-AdenoCa::TCGA-CD-A4MG-01A-11D-A25D-08",
  "Stomach-AdenoCa::TCGA-CD-A4MI-01A-21D-A25D-08",
  "Stomach-AdenoCa::TCGA-CD-A4MJ-01A-11D-A25D-08",
  "Stomach-AdenoCa::TCGA-CG-4305-01A-01D-1158-08",
  "Stomach-AdenoCa::TCGA-CG-4306-01A-01D-1158-08",
  "Stomach-AdenoCa::TCGA-CG-5721-01A-11D-1600-08",
  "Stomach-AdenoCa::TCGA-CG-5723-01A-11D-1600-08",
  "Stomach-AdenoCa::TCGA-CG-5726-01A-11D-1600-08",
  "Stomach-AdenoCa::TCGA-CG-5728-01A-11D-1600-08",
  "Stomach-AdenoCa::TCGA-CG-5733-01A-11D-1600-08",
  "Stomach-AdenoCa::TCGA-D7-A4YV-01A-11D-A25D-08",
  "Stomach-AdenoCa::TCGA-D7-A4YY-01A-11D-A25D-08",
  "Stomach-AdenoCa::TCGA-F1-6177-01A-11D-1800-08",
  "Stomach-AdenoCa::TCGA-F1-6874-01A-11D-1882-08",
  "Stomach-AdenoCa::TCGA-FP-A4BE-01A-12D-A24D-08",
  "Stomach-AdenoCa::TCGA-HF-7132-01A-11D-2053-08",
  "Stomach-AdenoCa::TCGA-HJ-7597-01A-21D-2201-08",
  "Stomach-AdenoCa::TCGA-HU-A4G8-01A-11D-A25D-08",
  "Stomach-AdenoCa::TCGA-HU-A4GN-01A-11D-A25D-08",
  "Stomach-AdenoCa::TCGA-HU-A4GQ-01A-11D-A25D-08",
  "Stomach-AdenoCa::TCGA-HU-A4GU-01A-11D-A25D-08",
  "Stomach-AdenoCa::TCGA-HU-A4GX-01A-12D-A25D-08",
  "Stomach-AdenoCa::TCGA-HU-A4H3-01A-21D-A25D-08",
  "Stomach-AdenoCa::TCGA-HU-A4H8-01A-11D-A25D-08",
  "Stomach-AdenoCa::TCGA-MX-A5UJ-01A-11D-A31L-08",
  "Stomach-AdenoCa::TCGA-R5-A7ZI-01A-11D-A34U-08",
  "Stomach-AdenoCa::TCGA-RD-A8NB-01A-12D-A397-08",
  "Stomach-AdenoCa::TCGA-SW-A7EA-01A-12D-A34U-08",
  "Stomach-AdenoCa::TCGA-VQ-A8P2-01A-11D-A364-08",
  "Stomach-AdenoCa::TCGA-VQ-A8PB-01A-11D-A397-08",
  "Stomach-AdenoCa::TCGA-VQ-A8PT-01A-11D-A410-08",
  "Stomach-AdenoCa::TCGA-VQ-A91D-01A-11D-A410-08",
  "Stomach-AdenoCa::TCGA-VQ-A91K-01A-11D-A410-08",
  "Stomach-AdenoCa::TCGA-VQ-A924-01A-11D-A410-08",
  "Stomach-AdenoCa::TCGA-VQ-AA6D-01A-11D-A410-08"
)

get_sig_universes <- function(exposures, filename, sigs, sig_file) {
  exposure_ids <- colnames(exposures)
  sample_types <- gsub("::.*$", "", exposure_ids)
  exposure_list <-
    PCAWG7::SplitMatrixBySampleType(exposures, sample_types)
  cancer_type_2_sig_prop <-
    lapply(exposure_list, function(exposure) {
      sig_prop <- apply(X = exposure, MARGIN = 1, FUN = function(counts) {
        length(counts[counts > 0]) / length(counts)
      })
      return(sig_prop[sig_prop > 0])
    })
  cat("", file = filename, append = FALSE)
  lapply(
    names(cancer_type_2_sig_prop),
    function(xx) {
      cat(xx, names(cancer_type_2_sig_prop[[xx]]),
        sep = ",", append = TRUE, file = filename
      )
      cat("\n", file = filename, append = TRUE)
      cat(xx, cancer_type_2_sig_prop[[xx]],
        sep = ",", append = TRUE, file = filename
      )
      cat("\n", file = filename, append = TRUE)
    }
  )
  sigs_to_use <- sigs[, rownames(exposures), drop = FALSE]
  ICAMS::WriteCatalog(
    catalog = sigs_to_use,
    file = sig_file
  )
  invisible(NULL)
}

get_sbs192_sig <- function() {
  tmp <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS192

  # Find out signatures which were previously extracted from exome data
  indices <- grep(pattern = "-E", x = colnames(tmp))

  to_keep <- tmp[, -indices]

  exome_sig <- tmp[, indices]

  # Change the abundance of exome signatures
  attributes(exome_sig)$region <- "exome"
  attributes(exome_sig)$abundance <-
    ICAMS::all.abundance$BSgenome.Hsapiens.1000genomes.hs37d5$exome$`192`

  # Transform the exome signatures to genome signatures
  genome_sig <-
    ICAMS::TransformCatalog(
      catalog = exome_sig,
      target.ref.genome = "GRCh37",
      target.region = "transcript",
      target.catalog.type = "counts.signature"
    )
  colnames(genome_sig) <-
    gsub(pattern = "-E", replacement = "", x = colnames(genome_sig))

  new_sig <- cbind(to_keep, genome_sig)
  new_sig2 <- new_sig[, gtools::mixedsort(colnames(new_sig))]
  return(new_sig2)
}

get_pcawg_sbs192_exposure <- function() {
  sbs96_spectra <- PCAWG7::spectra$PCAWG$SBS96
  sample_names <- colnames(sbs96_spectra)

  sbs192_spectra <- PCAWG7::spectra$PCAWG$SBS192[, sample_names]
  sbs96_exposure <- PCAWG7::exposure$PCAWG$SBS96[, sample_names]

  sbs96_mutation <- colSums(sbs96_spectra)
  sbs192_mutation <- colSums(sbs192_spectra)

  scale_factor <- sbs192_mutation / sbs96_mutation

  retval <- lapply(sample_names, FUN = function(sample_name) {
    sample_sbs96_exposure <- sbs96_exposure[, sample_name, drop = FALSE]
    sample_scale_factor <- scale_factor[sample_name]
    sample_sbs192_exposure <- sample_sbs96_exposure * sample_scale_factor
    return(round(sample_sbs192_exposure))
  })

  sbs192_exposure <- do.call("cbind", retval)

  return(sbs192_exposure)
}

check_syn_exposures <- function(exposure) {
  msi_sig_names <-
    c("SBS6", "SBS14", "SBS15", "SBS20", "SBS21", "SBS26", "SBS44")
  tmp <- 
    exposure[rownames(exposure) %in% msi_sig_names, , drop = FALSE]
  tmp2 <- tmp[, colSums(tmp) > 0, drop = FALSE]
  names_to_change <- colnames(tmp2)
  indices <- which(colnames(exposure) %in% names_to_change)
  colnames(exposure)[indices] <-
    gsub(pattern = "S.", replacement = "MSI-H.S.", 
         x = colnames(exposure)[indices])
  return(exposure)
}