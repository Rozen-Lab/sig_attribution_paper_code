library(ICAMS)
library(data.table)

get_all_input <- function(dataset_name,
                          data_top_folder_name = "synthetic_data") {
  data_home <- file.path(data_top_folder_name, dataset_name)
  
  spectra <- ICAMS::ReadCatalog(
    file = file.path(data_home, "ground.truth.syn.catalog.csv")
  )
  
  ground_truth_sigs <- ICAMS::ReadCatalog(
    file.path(data_home, "ground.truth.sigs.csv"),
    catalog.type = "counts.signature"
  )
  
  sig_universe_info <- data.table::fread(
    file.path(data_home, "ground.truth.sig.universe.csv"),
    fill = TRUE
  )
  
  spectra_list <- PCAWG7::SplitPCAWGMatrixByTumorType(spectra)
  
  cancer_types <- names(spectra_list)
  
  signature_universes <-
    sapply(cancer_types, FUN = function(cancer_type) {
      one_type_info <- sig_universe_info[V1 == cancer_type, -1]
      sig_props <- unlist(one_type_info[2, ])
      sig_props <- as.numeric(sig_props[sig_props != ""])
      
      sig_names <- unlist(one_type_info[1, ])
      sig_names <- sig_names[sig_names != ""]
      names(sig_props) <- sig_names
      return(sig_props)
    })
  
  return(list(
    spectra_list = spectra_list,
    ground_truth_sigs = ground_truth_sigs,
    signature_universes = signature_universes,
    cancer_types = cancer_types,
    all_spectra  = spectra
  ))
}

get_ground_truth_exposure = function(mutation_type,
                                     data_top_folder_name = "synthetic_data") {
  mSigTools::read_exposure(
    file.path(data_top_folder_name, 
              mutation_type, 
              "ground.truth.syn.exposures.csv"))
}
