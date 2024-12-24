
source("analysis/code/get_all_input.R")
library(cosmicsig)
library(parallel)

reoptimize_similarities_from_exposures =
  function(
    spectra, 
    exposures, 
    signatures = 
      cbind(
        cosmicsig::COSMIC_v3.4$signature$GRCh37$SBS96,
        cosmicsig::COSMIC_v3.3$signature$GRCh37$SBS96[ , c("SBS22", "SBS40"), drop = FALSE])
  )
  {
    stopifnot(ncol(spectra) == ncol(exposures))
    
    do1index = function(ii) {
      message(ii)
      retval = reoptimize_exposure_w_ll(
        spectra[ , ii, drop = FALSE],
        exposures[ , ii, drop = FALSE],
        signatures)
    }
    
    retlist = mclapply(1:ncol(spectra), 
                       do1index,
                       mc.cores = 100)
    return(mSigAct:::ListOfList2Tibble(retlist))
    
  }

reoptimize_exposure = function(target_spectrum, source_exposure, sig_matrix) {
  # browser()
  usedsigs <- rownames(source_exposure)[source_exposure[, 1] > 0]
  missing_sigs = setdiff(usedsigs, colnames(sig_matrix))
  if (length(missing_sigs) != 0) {
    stop(paste("signature(s)", paste(missing_sigs, collapse = " "),
               "are missing in sample", colnames(target_spectrum)))
  }
  
  reduced_sig_matrix <- as.matrix(sig_matrix)[ , usedsigs, drop = FALSE]
  # message("ncol(nnls reduced_sig_matrix = ", ncol(reduced_sig_matrix))
  reoptimized_exposure <- 
    mSigTools::optimize_exposure_QP(target_spectrum, reduced_sig_matrix)
  recon <- reduced_sig_matrix %*% reoptimized_exposure
  recon = round(recon)
  reoptimized_sim <- mSigAct::cossim(recon, target_spectrum)

  L2 = philentropy::distance(
    t(cbind(recon,
            target_spectrum[, 1])),
    method = 'squared_euclidean')
  L2 = L2 / sum(target_spectrum)
  md = philentropy::distance(
    t(cbind(recon,
            target_spectrum[ , 1])),
    method = 'manhattan')
  md = md/sum(target_spectrum)

  return(list(
    sample_id            = colnames(target_spectrum),
    cossim_for_NNLS      = reoptimized_sim,
    exposure_for_NNLS    = reoptimized_exposure,
    reduced_sig_matrix   = reduced_sig_matrix,
    smd_for_NNLS = md,
    L2_for_NNLS = L2
  ))
}


reoptimize_exposure_w_ll =
  function(target_spectrum, source_exposure, sig_matrix) {
    list1 = reoptimize_exposure(
      target_spectrum, 
      source_exposure, 
      sig_matrix)
    # browser()
    reduced_sig_matrix = list1$reduced_sig_matrix
    llopts = mSigAct::DefaultManyOpts('multinom')
    # message("ncol(ll reduced_sig_matrix = ", ncol(reduced_sig_matrix))
    retval2 = mSigAct:::OptimizeExposure(
      spectrum = target_spectrum,
      sigs     = reduced_sig_matrix,
      m.opts   = llopts)
    # browser()
    optimized_exposure = retval2$exposure
    recon = reduced_sig_matrix %*% optimized_exposure
    recon = round(recon)
    L2 = philentropy::distance(
      t(cbind(recon,
              target_spectrum[ , 1])),
      method = 'squared_euclidean')
    L2 = L2 / sum(target_spectrum)
    md = philentropy::distance(
      t(cbind(recon,
              target_spectrum[ , 1])),
      method = 'manhattan')
    md = md / sum(target_spectrum)
    # browser()
    reoptimized_sim = 
      mSigAct::cossim(recon, target_spectrum)
    list1$reduced_sig_matrix = NULL
    return(
      c(
        list1,
        list(logLH = retval2$loglh,
             exposure_for_LH = optimized_exposure,
             cossim_for_LH = reoptimized_sim,
             smd_for_LH = md,
             L2_for_LH = L2
        )))
  }


pcawg_stomach_spectra = function() {
  s96 <- PCAWG7::spectra$PCAWG$SBS96
  all_stomach_spectra <-
    s96[, grep("Stomach-AdenoCA", colnames(s96), fixed = TRUE)]
  return(all_stomach_spectra)
}


pcawg_stomach_exposures = function() {
  real_exp = PCAWG7::exposure$PCAWG$SBS96
  real_exp = real_exp[ , grep("Stomach-AdenoCA", colnames(real_exp))]
  stopifnot(dim(real_exp) == c(65, 75))
  stopifnot(sum(rowSums(real_exp) > 0) == 20)
  return(real_exp)
}


synthetic_stomach_exposures = function() {
  real_exp = get_ground_truth_exposure("SBS")
  real_exp = real_exp[ , grep("Stomach-AdenoCA", colnames(real_exp))]
  # Not all exposures are present in stomach, however
  stopifnot(dim(real_exp) == c(41, 100))
  return(real_exp)
}


synthetic_stomach_spectra = function() {
  ret1 = get_all_input("SBS")
  # browser()
  all_spectra = ret1$all_spectra
  stomach_spectra = 
    all_spectra[ , 
                 grep("Stomach-AdenoCA", 
                      colnames(all_spectra))]
  stopifnot(dim(stomach_spectra) == c(96, 100))
  return(stomach_spectra)
}
