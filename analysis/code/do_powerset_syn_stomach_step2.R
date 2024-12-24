if(!grepl("attribution_paper_code", basename(getwd()))) {
  stop("Run this script from the top level directory")
}
rm(list = ls())
message("ls() = ", ls())
source("analysis/code/powerset_analysis.R")
library(huxtable)
library(openxlsx)


## OLD all_stomach_spectra <- pcawg_stomach_spectra()

# if no git LSF
# bigfile = "c:/users/steve rozen/Documents/bigdata/syn_stomach_powerset.rdata"
# if git LSF 
bigfile = "output_for_paper/syn_stomach_powerset.rdata"

stopifnot(file.exists(bigfile))
message("Loading existing file ", bigfile)
load(bigfile)
stopifnot(length(syn_stomach_powerset) == 100)

# Compare ground truth attribution to attributions in the power set
# Used for checking how many attributions have higher cosine
# similarity than the ground truth attributions
compare_to_powerset <- function(powerset_info, syn_reoptimized, exposures, 
                                sim_threshold = 0.969,
                                proportion_threshold = 0.03) {

  # browser()
  sampleid <- powerset_info[[1]]
  rest <- powerset_info[[2]]
  # target_spectrum <- spectra[, sampleid, drop = FALSE]
  ground_truth_exposure <- exposures[ , sampleid, drop = FALSE]
  

  usedsigs <- rownames(ground_truth_exposure)[ground_truth_exposure[, 1] > 0]

  gt_exposure_used = 
    unlist(syn_reoptimized[syn_reoptimized$sample_id == sampleid, ]$exposure_for_NNLS)
    
  reoptimized_sim = 
    dplyr::filter(syn_reoptimized, sample_id == sampleid)$cossim_for_NNLS
  
  ground_truth_num_sigs <- length(usedsigs)

  processexps <- function(zz) {
    cosine <- zz[[4]]
    return(cosine > reoptimized_sim)
  }
  
  processexps_filtered_by_proportion <- function(zz) {
    cosine <- zz[[4]]
    if (!filter_exposure_no_low_fraction_contribution(zz, proportion_threshold)) {
      return(FALSE)
    }
    return(cosine > reoptimized_sim)
  }
  
  cosines_filtered_by_proportion <- function(zz) {
    if (!filter_exposure_no_low_fraction_contribution(zz, proportion_threshold)) {
      return(0)
    }
    return(zz[[4]])
  }
  
  exps <- rest[["exposures"]]
  # exps is a list of the powerset exposures that generated
  # reconstructions that passed the cosine similarity 
  # threshold.  The exposure is in slot q_exp_matrix of
  # each element of exp. The similarity of the reconstruction
  # is in slot 4, which doesn't have a name
  
  num_better =
    sum(unlist(lapply(exps, processexps)))

  tmp = unlist(lapply(exps, processexps_filtered_by_proportion))
  
  filtered_better_indices = which(tmp)
  
  filtered_num_better = length(filtered_better_indices)
  # browser()
  ok_cosines = unlist(lapply(exps, cosines_filtered_by_proportion))

  best_ok_cosine = max(ok_cosines)
  best = which(ok_cosines == best_ok_cosine)
  expinfo = exps[[best[1]]]
  bestsigs = rownames(expinfo$q_exp_matrix)
  
  ground_truth_exp_ok = 
    exposure_exceeds_proportion_filter(
      gt_exposure_used,
      proportion_threshold = proportion_threshold)
  
  return(list(
    sample_id = sampleid,
    ground_truth_sigs = usedsigs,
    ground_truth_sigs_string = paste(usedsigs, collapse = " "),
    ground_truth_similarity = reoptimized_sim,
    ground_truth_similarity_ok = reoptimized_sim > sim_threshold,
    ground_truth_exposure_ok = ground_truth_exp_ok,
    number_better_attributions = num_better,
    number_better_attributions_filtered_by_proportion = filtered_num_better,
    best_ok_alternative_similarity = best_ok_cosine,
    best_attribution = bestsigs,
    best_attribution_sigs = paste(bestsigs, collapse = " ")
  ))
}

get_alter_attribution_table <- function(powerset_results, 
                                        syn_reoptimized,
                                        gt_exposures) {
  all.morexx <- lapply(powerset_results, 
                       compare_to_powerset, 
                       syn_reoptimized = syn_reoptimized,
                       exposures = gt_exposures)

  mSigAct:::ListOfList2Tibble(all.morexx) %>%
    dplyr::mutate(sample_id = gsub("Stomach-AdenoCA::", "", sample_id,
                                   fixed = TRUE)) ->
    table1
  
  more = apply(table1, 1, function(x) {
    FN = setdiff(unlist(x$ground_truth_sigs), unlist(x$best_attribution))
    FP = setdiff(unlist(x$best_attribution), unlist(x$ground_truth_sigs))
    list(
      FN_string = paste(FN, collapse = " "),
      num_FN = length(FN),
      FP_string = paste(FP, collapse = " "),
      num_FP = length(FP),
      correct = (length(FN) == 0 & length(FP) == 0)
    )
  })
  more2 = mSigAct:::ListOfList2Tibble(more)
  
  table2 = cbind(table1, more2)
 
  dplyr::mutate(
    table2,
      best_attribution = NULL,
    ground_truth_sigs = ground_truth_sigs_string) ->
    table3
  
  return(table3)
}




table_results <- 
  get_alter_attribution_table(syn_stomach_powerset,
                              syn_all_reoptimized,
                              get_ground_truth_exposure("SBS"))

sum(table_results$correct)
# 13
sum(!table_results$correct & table_results$ground_truth_exposure_ok)
# 17
sum(!table_results$correct & !table_results$ground_truth_exposure_ok)
# 70
mean(table_results$num_FP)
# 0.82
mean(table_results$num_FN)
# [1] 1.17


sum(table_results$number_better_attributions > 0) # 99
median(table_results$number_better_attributions) # 37

table_results %>%
  dplyr::mutate(ground_truth_sigs_string = NULL,
                FN_string = NULL,
                FP_string = NULL,
                ground_truth_similarity_ok = NULL
                ) %>%
  huxtable::as_huxtable(add_rownames = FALSE, add_colnames = TRUE) %>%
  insert_row(
    paste("Table S5: For Stomach-AdenoCA synthetic data,",
          "the cosine similarity",
          "XXXX",
          "XXXX"),
    fill = "",
    after = 0) %>%
  set_wrap(row = 1, value = TRUE) %>%
  merge_cells(1, 1:2) %>%
  set_header_rows(1, TRUE) %>%
  style_header_rows(bold = TRUE) %>%
  set_align(2:nrow(.), 1:ncol(.), "center") %>%
  set_number_format(2:nrow(.), c(3,7), "%5.4f") %>%
  huxtable::quick_xlsx(
    file = "output_for_paper/table_s5_synthetic_stomach_powerset.xlsx")
