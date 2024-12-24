convert_to_msa_format <- function(matrix, mut_type) {
  df <- as.data.frame(matrix)
  old_row_names <- rownames(df)
  if (mut_type %in% c("DBS", "DBS78")) {
    msa_row_names <- paste0(
      substr(old_row_names, 1, 2), ">",
      substr(old_row_names, 3, 4)
    )
  } else if (mut_type == "ID") {
    msa_row_names <-
      gsub(pattern = ":", replacement = "_", x = old_row_names)
  }

  msa_df <- cbind("Mutation type" = msa_row_names, df)
  return(msa_df)
}

write_msa_file <-
  function(catalog, sig, output_dir, mut_type, data_identifier) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    if (mut_type %in% c("SBS", "SBS96")) {
      catalog_file_name <- paste0("WGS_", data_identifier, ".96.csv")
      ICAMS::WriteCatalog(
        catalog = catalog,
        file = file.path(output_dir, catalog_file_name)
      )

      sig_file_name <- paste0(data_identifier, "_SBS_signatures.csv")
      ICAMS::WriteCatalog(
        catalog = sig,
        file = file.path(output_dir, sig_file_name)
      )
    } else if (mut_type %in% c("DBS", "DBS78")) {
      msa_catalog <-
        convert_to_msa_format(matrix = catalog, mut_type = mut_type)
      catalog_file_name <- paste0("WGS_", data_identifier, ".dinucs.csv")
      write.csv(
        x = msa_catalog,
        file = file.path(output_dir, catalog_file_name), row.names = FALSE
      )

      msa_sig <-
        convert_to_msa_format(matrix = sig, mut_type = mut_type)
      sig_file_name <- paste0(data_identifier, "_DBS_signatures.csv")
      write.csv(
        x = msa_sig,
        file = file.path(output_dir, sig_file_name), row.names = FALSE
      )
    } else if (mut_type == "ID") {
      msa_catalog <-
        convert_to_msa_format(matrix = catalog, mut_type = mut_type)
      catalog_file_name <- paste0("WGS_", data_identifier, ".indels.csv")
      write.csv(
        x = msa_catalog,
        file = file.path(output_dir, catalog_file_name), row.names = FALSE
      )

      msa_sig <-
        convert_to_msa_format(matrix = sig, mut_type = mut_type)
      sig_file_name <- paste0(data_identifier, "_ID_signatures.csv")
      write.csv(
        x = msa_sig,
        file = file.path(output_dir, sig_file_name), row.names = FALSE
      )
    }
  }

write_msi_sample_sig_for_msa <-
  function(spectra, msi_samples, sig, output_dir, mut_type, data_identifier) {
    output_dir_non_msi <- file.path(output_dir, "non_msi")
    output_dir_msi <- file.path(output_dir, "msi")

    if (length(msi_samples) == 0) {
      write_msa_file(
        catalog = spectra, sig = sig,
        output_dir = file.path(output_dir_non_msi, data_identifier),
        mut_type = mut_type, data_identifier = data_identifier
      )
    } else {
      spectra_msi <- spectra[, msi_samples, drop = FALSE]
      sig_msi <- sig
      write_msa_file(
        catalog = spectra_msi, sig = sig_msi,
        output_dir = file.path(output_dir_msi, data_identifier),
        mut_type = mut_type, data_identifier = data_identifier
      )

      non_msi_samples <- setdiff(colnames(spectra), msi_samples)
      spectra_non_msi <- spectra[, non_msi_samples, drop = FALSE]
      msi_sig_names <- get_msi_sig_names()
      non_msi_sig_names <- setdiff(colnames(sig), msi_sig_names)
      sig_non_msi <- sig[, non_msi_sig_names, drop = FALSE]

      write_msa_file(
        catalog = spectra_non_msi, sig = sig_non_msi,
        output_dir = file.path(output_dir_non_msi, data_identifier),
        mut_type = mut_type, data_identifier = data_identifier
      )
    }
  }

classify_msi_sample_sig_syn_for_msa <-
  function(dataset_name, output_home, data_identifier,
           data_top_folder_name = "synthetic_data",
           cancer_types = NULL) {
    all_inputs <- get_all_input(
      dataset_name = dataset_name,
      data_top_folder_name = data_top_folder_name
    )
    all_spectra <- all_inputs$spectra_list
    gt_sig <- all_inputs$ground_truth_sigs
    sig_universe <- all_inputs$signature_universes

    if (is.null(cancer_types)) {
      cancer_types <- all_inputs$cancer_types
    }

    for (cancer_type in cancer_types) {
      spectra <- all_spectra[[cancer_type]]
      sig_names <- names(sig_universe[[cancer_type]])
      sig <- gt_sig[, sig_names, drop = FALSE]
      output_dir <- file.path(output_home, cancer_type)

      # Check whether there are any MSI-H samples in spectra
      msi_samples <- grep(pattern = "MSI", x = colnames(spectra), value = TRUE)
      write_msi_sample_sig_for_msa(
        spectra = spectra, msi_samples = msi_samples,
        sig = sig, output_dir = output_dir,
        mut_type = dataset_name, data_identifier = data_identifier
      )
    }
  }

prepare_msa_bash_script <-
  function(bash_script_dir, bash_file_name, msa_nextflow_file,
           data_dir, data_name,
           output_dir, type,
           project_dir = getwd(), tool = "conda") {
    if (!dir.exists(bash_script_dir)) {
      dir.create(path = bash_script_dir, recursive = TRUE)
    }
    if (type == "SBS96") {
      type <- "SBS"
    }

    if (type == "DBS78") {
      type <- "DBS"
    }

    file_name <- file.path(bash_script_dir, paste0(bash_file_name, ".txt"))
    my_file <- file(description = file_name, open = "w")

    writeLines("#!/bin/bash", my_file)
    writeLines(paste0("PROJECT_DIR=", project_dir), my_file)
    writeLines(paste0("MSA_NF_FILE=", msa_nextflow_file), my_file)
    writeLines(paste0("DATA_DIR=", data_dir), my_file)
    writeLines(paste0("DATA_NAME=", data_name), my_file)
    writeLines(paste0("OUTPUT_DIR=", output_dir), my_file)
    writeLines(paste0("TOOL=", tool), my_file)
    writeLines(paste0("TYPE=", type), my_file)
    writeLines("", my_file)
    writeLines('ABS_MSA_NF_FILE="$PROJECT_DIR"/$MSA_NF_FILE', my_file)
    writeLines('ABS_DATA_DIR="$PROJECT_DIR"/$DATA_DIR', my_file)
    writeLines('ABS_OUTPUT_DIR="$PROJECT_DIR"/$OUTPUT_DIR', my_file)
    writeLines("", my_file)
    writeLines("nice nextflow run $ABS_MSA_NF_FILE \\", my_file)
    writeLines(" -profile $TOOL --dataset $DATA_NAME --input_tables $ABS_DATA_DIR \\", my_file)
    writeLines(" --mutation_types $TYPE --output_path $ABS_OUTPUT_DIR \\", my_file)
    writeLines(" --signature_prefix $DATA_NAME --signature_tables $ABS_DATA_DIR/$DATA_NAME \\", my_file)
    writeLines("", my_file)
    writeLines('TAR_FILE_DIR="$ABS_OUTPUT_DIR"/output_tables/', my_file)
    writeLines('TAR_FILE_NAME=MSA_output_"$DATA_NAME".tar.gz', my_file)
    writeLines("cd $TAR_FILE_DIR", my_file)
    writeLines("tar xf $TAR_FILE_NAME", my_file)
    writeLines('EXPOSURE_FILE_DIR="$TAR_FILE_DIR"/$DATA_NAME', my_file)
    writeLines('EXPOSURE_FILE_NAME=output_"$DATA_NAME"_"$TYPE"_mutations_table.csv', my_file)
    writeLines('PRUNED_EXPOSURE_FILE_NAME=pruned_attribution_"$DATA_NAME"_"$TYPE"_abs_mutations.csv', my_file)
    writeLines("cp $EXPOSURE_FILE_DIR/$EXPOSURE_FILE_NAME $ABS_OUTPUT_DIR/..", my_file)
    writeLines("cp $EXPOSURE_FILE_DIR/$PRUNED_EXPOSURE_FILE_NAME $ABS_OUTPUT_DIR/..", my_file)
    writeLines("cp $ABS_OUTPUT_DIR/nf-pipeline_info/MSA-nf_report.html $ABS_OUTPUT_DIR/..", my_file)

    close(my_file)
    files <- list.files(
      path = bash_script_dir,
      pattern = ".txt", full.names = TRUE
    )
    newfiles <- gsub(".txt$", ".sh", files)
    file.rename(files, newfiles)
  }

prepare_msa_data_script_syn <-
  function(dataset_name, input_dir, data_identifier,
           output_home, bash_script_home, msa_nextflow_file,
           data_top_folder_name = "synthetic_data",
           cancer_types = NULL,
           project_dir = getwd()) {
    classify_msi_sample_sig_syn_for_msa(
      dataset_name = dataset_name,
      output_home = input_dir,
      data_identifier = data_identifier,
      data_top_folder_name = data_top_folder_name,
      cancer_types = cancer_types
    )

    cancer_types <-
      list.dirs(input_dir, recursive = FALSE, full.names = FALSE)

    bash_script_dir <-
      file.path(bash_script_home, "bash", dataset_name, "syn")

    for (cancer_type in cancer_types) {
      sub_folders <- list.files(file.path(input_dir, cancer_type))

      for (sub_folder in sub_folders) {
        full_input_dir <- file.path(input_dir, cancer_type, sub_folder)
        output_dir <- file.path(output_home, cancer_type, sub_folder, "raw")
        bash_file_name <-
          paste0("run_", cancer_type, "_", sub_folder, "_", dataset_name)
        prepare_msa_bash_script(
          bash_script_dir = bash_script_dir,
          bash_file_name = bash_file_name,
          msa_nextflow_file = msa_nextflow_file,
          data_dir = full_input_dir,
          data_name = data_identifier,
          output_dir = output_dir,
          type = dataset_name,
          project_dir = project_dir
        )
      }
    }

    bash_files <- list.files(path = bash_script_dir, pattern = ".sh")

    file_name <-
      file.path(bash_script_home, paste0("run_bash_script_syn_", dataset_name, ".txt"))
    my_file <- file(description = file_name, open = "w")

    writeLines("#!/bin/bash", my_file)
    writeLines(paste0("PROJECT_DIR=", project_dir), my_file)
    writeLines(paste0("BASH_SCRIPT_DIR=", bash_script_dir), my_file)
    writeLines('ABS_BASH_SCRIPT_DIR="$PROJECT_DIR"/$BASH_SCRIPT_DIR', my_file)
    writeLines("cd $ABS_BASH_SCRIPT_DIR", my_file)
    writeLines("", my_file)
    for (bash_file in bash_files) {
      writeLines(paste0("bash ", bash_file), my_file)
    }
    close(my_file)
    files <- list.files(
      path = bash_script_home,
      pattern = ".txt", full.names = TRUE
    )
    newfiles <- gsub(".txt$", ".sh", files)
    file.rename(files, newfiles)
  }



summarize_msa_results <- function(msa_output_dir, mut_type) {
  pruned_exposure_files <-
    list.files(
      path = msa_output_dir, pattern = "pruned_attribution",
      full.names = TRUE, recursive = TRUE
    )

  pruned_exposure_files2 <-
    grep(
      pattern = "/raw/", x = pruned_exposure_files,
      invert = TRUE, value = TRUE
    )

  pruned_msa_exposures <-
    lapply(pruned_exposure_files2, FUN = function(file) {
      t(read.table(file, header = TRUE, sep = ",", row.names = 1))
    })

  pruned_msa_exposure_all <-
    mSigAct:::MergeListOfExposures(pruned_msa_exposures)

  mSigTools::write_exposure(
    exposure = pruned_msa_exposure_all,
    file = file.path(msa_output_dir, "inferred_exposures.csv")
  )
}
