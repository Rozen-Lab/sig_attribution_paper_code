#!/usr/bin/env nextflow

// Copyright (C) 2022 Sergey Senkin

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// input data
params.SP_extractor_output_path = null // optional path to SigProfilerExtractor output
params.SP_matrix_generator_output_path = null // optional path to SigProfilerMatrixGenerator output
params.COSMIC_signatures = false // if set to true, COSMIC signatures are used form SigProfiler output, otherwise de-novo ones are used
params.dataset = ['SIM_test'] // several datasets can be provided as long as input mutation tables are available
params.mutation_types = ['SBS', 'DBS', 'ID'] // add or remove mutation types if needed
params.input_tables = "$baseDir/input_mutation_tables"
params.SBS_context = 96 // 96, 192, 288, 1536 context matrices can be provided (SBS only)
params.number_of_samples = -1 // number of samples to analyse (-1 means all available)

// output paths
params.output_path = "."
params.tables_output_path = params.output_path + "/output_tables"
params.plots_output_path = params.output_path + "/plots"

// signatures to use
params.signature_tables = "$baseDir/signature_tables"
params.signature_prefix = "sigProfiler" // prefix of signature files to use (e.g. sigProfiler, sigRandom)

// simulations parameters
params.run_only_simulations = false // set to true if only simulations are needed, these will be produced in $baseDir/output_tables folder
params.number_of_simulated_samples = -1 // number of simulations to run (-1 means automatically apply a rounded factor of ten but not less than 1000)
params.add_noise = true // add noise in simulations (recommended)
params.noise_type = "gaussian" // set the type of noise in simulations: gaussian, poisson or negative_binomial (Gaussian by default)
params.noise_stdev = 10 // set standard deviation of gaussian noise, in percentage of sample mutation burden (10 percent by default)
params.zero_inflation_threshold = 0.05 // set the relative threshold below which all simulated signature activities are set to zero (0.01 by default)

// optimisation flag and parameters
params.run_only_optimisation = false // set to true if only optimisation is required, without final attributions
params.optimisation_NNLS_output_path = "$baseDir/outputs_optimisation"
params.optimisation_plots_output_path = params.plots_output_path + "/optimisation_plots"
params.optimised = true // if set to false, optimisation will run but not be used in final attributions
params.optimisation_strategy = "removal" // optimisation strategy (removal, addition or add-remove)
params.weak_thresholds = ['0.0000', '0.0001', '0.0002', '0.0003'] // range of L2 similarity decrease thresholds to be scanned, excluding weakest signatures - adjust if needed
params.strong_thresholds = ['0.0000'] // range of L2 similarity increase thresholds to be scanned, including strongest signatures: only one is sufficient in default removal strategy
params.bootstrap_method = "binomial" // bootstrap flag and method (binomial, multinomial, residuals, classic, bootstrap_residuals)
params.number_of_bootstrapped_samples_in_optimisation = 100 // at least 100 is recommended
params.number_of_bootstrapped_samples = 1000 // bootstrap variations in final attribution, at least 1000 is recommended
params.confidence_level = 0.95 // specify the confidence level for CI calculation (default: 0.95)
params.use_absolute_attributions = false // use absolute mutation counts in final bootstrap outputs (relative by default)
params.metric_to_prioritise = "specificity" // set a metric to prioritise (default: specificity), requiring at least the specified threshold or closest alternative
params.metric_threshold = 0.95 // specify the minimum threshold of the prioritised metric
params.signatures_to_prioritise = [] // set a list of signatures to prioritise (empty list means all, by default)
params.no_CI_for_penalties = false // do not use confidence intervals for optimal penalties calculation
params.calculate_penalty_on_average = false // apply criteria based on signatures overall (on average, less conservative), rather than maximising prioritised metric for every signature (more conservative)

// plotting flags
params.plot_optimisation_plots = true
params.plot_bootstrap_attributions = true
params.plot_metrics = true
params.plot_signatures = true
params.plot_input_spectra = true
params.plot_fitted_spectra = false
params.plot_residuals = false
params.show_poisson_errors = false
params.show_strands = false // only works with higher contexts (192, 288)
params.show_nontranscribed_region = false // only wortks with higher contexts (288)

// if SIM in dataset name (synthetic data), use the following percentage range for measuring signature attirbution sensitivities
params.signature_attribution_thresholds = 0..20

// helper flags for scripts (automatic based on parameters)
optimised_flag = (params.optimised) ? "-x" : ''
abs_flag = (params.use_absolute_attributions) ? "-a" : ''
suffix = (params.use_absolute_attributions) ? "abs_mutations" : 'weights'
error_flag = (params.show_poisson_errors) ? "-e" : ''
strands_flag = (params.show_strands) ? "-b" : ''
nontranscribed_flag = (params.show_nontranscribed_region) ? "-n" : ''
COSMIC_flag = (params.COSMIC_signatures) ? "-C" : ''
signature_prefix = (params.SP_extractor_output_path) ? params.signature_prefix + "_conv" : params.signature_prefix
noise_flag = (params.add_noise) ? "-z" : ''
no_CI_for_penalties_flag = (params.no_CI_for_penalties) ? "--no_CI" : ''
calculate_penalty_on_average_flag = (params.calculate_penalty_on_average) ? "--average" : ''
prioritised_signatures_flag = (params.signatures_to_prioritise) ? "--signatures_to_prioritise " + params.signatures_to_prioritise.join(' ') : ''
// override number of samples/variations for test run
test_run = ((params.dataset == ['SIM_test']) || (params.dataset == 'SIM_test')) ? true : false
number_of_bootstrapped_samples_in_optimisation = (test_run) ? 10 : params.number_of_bootstrapped_samples_in_optimisation
number_of_bootstrapped_samples = (test_run) ? 10 : params.number_of_bootstrapped_samples
number_of_simulated_samples = (test_run) ? 10 : params.number_of_simulated_samples
weak_thresholds = (test_run) ? ['0.0000', '0.0100', '0.0200'] : params.weak_thresholds
strong_thresholds = (test_run) ? ['0.0000'] : params.strong_thresholds

params.help = null

log.info ''
log.info '--------------------------------------------------------'
log.info '              __   __  _____                             '
log.info '             |  \\/  |/ ____|  /\\                      '
log.info '             | \\  / | (___   /  \\                     '
log.info '             | |\\/| |\\___ \\ / /\\ \\                 '
log.info '             | |  | |____) / ____ \\                    '
log.info '             |_|  |_|_____/_/    \\_\\                  '
log.info '                                                        '
log.info '          MUTATIONAL SIGNATURE ANALYSIS v2.2            '
log.info '--------------------------------------------------------'
log.info 'Copyright (C) Sergey Senkin'
log.info 'This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE'
log.info 'This is free software, and you are welcome to redistribute it'
log.info 'under certain conditions; see LICENSE for details.'
log.info '--------------------------------------------------------'
log.info ''

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run run_analysis.nf"
    log.info ""
    log.info "Nextflow currently does not support list parameters,"
    log.info "so please specify the parameters directly in the script."
    log.info ""
    exit 0
} else {
/* Software information */
log.info "help:                               ${params.help}"
}

// add parameter values to log output (.nextflow.log)
log.info params.collect { k,v -> "${k.padRight(34)}: $v" }.join("\n")

if (params.SP_extractor_output_path) {
  process convert_signature_tables {
    publishDir "${params.signature_tables}", mode: 'move', overwrite: true

    input:
    path input_path from params.SP_extractor_output_path

    output:
    file '*.csv' into signatures_for_spectra
    file '*.csv' into signatures_for_unoptimised_NNLS

    script:
    """
    python $baseDir/bin/convert_SP_to_MSA.py -S -t ${params.mutation_types.join(' ')} \
                                             -n ${signature_prefix} ${COSMIC_flag} \
                                             -i ${input_path} -s ${params.signature_tables} -o "./"
    """
  }
  if (!params.SP_matrix_generator_output_path) {
    process convert_extractor_input_matrices {
      publishDir "${params.input_tables}", mode: 'move', overwrite: true

      input:
      each dataset from params.dataset
      path input_path from params.SP_extractor_output_path

      output:
      file '*/*.csv' into converted_SP_to_MSA_for_spectra
      file '*/*.csv' into converted_SP_to_MSA_for_unoptimised_NNLS

      script:
      """
      python $baseDir/bin/convert_SP_to_MSA.py -E -d ${dataset} -t ${params.mutation_types.join(' ')} \
                                              -i ${input_path} -s ${params.signature_tables} -o "./"
      """
    }
  }
} else {
  // placeholder channels for execution from existing input matrices
  signatures_for_spectra = Channel.value(1)
  signatures_for_unoptimised_NNLS = Channel.value(1)
}

if (params.SP_matrix_generator_output_path) {
  process convert_input_data {
    publishDir "${params.input_tables}", mode: 'move', overwrite: true

    input:
    each dataset from params.dataset
    path input_path from params.SP_matrix_generator_output_path

    output:
    file '*/*.csv' into converted_SP_to_MSA_for_spectra
    file '*/*.csv' into converted_SP_to_MSA_for_unoptimised_NNLS

    script:
    """
    python $baseDir/bin/convert_SP_to_MSA.py -d ${dataset} -t ${params.mutation_types.join(' ')} \
                                             -i ${input_path} -s ${params.signature_tables} -o "./"
    """
  }
} else if (!params.SP_extractor_output_path) {
  // placeholder channels for execution from existing input matrices
  converted_SP_to_MSA_for_spectra = Channel.value(1)
  converted_SP_to_MSA_for_unoptimised_NNLS = Channel.value(1)
}

process plot_input_spectra {
  tag "${mutation_type}/${dataset}"
  publishDir "${params.plots_output_path}", mode: 'move'

  input:
  each mutation_type from params.mutation_types
  each dataset from params.dataset
  file inputs from converted_SP_to_MSA_for_spectra

  output:
  file '*/*/*.pdf' optional true
  file '*/*/*/*.pdf' optional true
  file '*/*/*/*/*.pdf' optional true
  file '*/*/*/*/*/*.pdf' optional true

  when:
  params.plot_input_spectra

  script:
  """
  python $baseDir/bin/plot_mutation_spectra.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                -i ${params.input_tables} ${strands_flag} ${nontranscribed_flag} -o "./"
  python $baseDir/bin/plot_mutation_spectra.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                -r -i ${params.input_tables} ${strands_flag} ${nontranscribed_flag} -o "./"
  """
}

process plot_signatures {
  tag "${mutation_type}/${dataset}"
  publishDir "${params.plots_output_path}", mode: 'move'

  input:
  each mutation_type from params.mutation_types
  each dataset from params.dataset
  file signatures from signatures_for_spectra

  output:
  file '*/*/*.pdf' optional true
  file '*/*/*/*.pdf' optional true
  file '*/*/*/*/*.pdf' optional true
  file '*/*/*/*/*/*.pdf' optional true

  when:
  params.plot_signatures

  script:
  """
  python $baseDir/bin/plot_mutation_spectra.py -S -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} \
                                               -p ${signature_prefix} -s ${params.signature_tables} \
                                               -r ${strands_flag} ${nontranscribed_flag} -o "./"
  """
}

process run_unoptimised_model {
  tag "${mutation_type}/${dataset}"
  publishDir "$baseDir/output_tables_unoptimised", mode: 'copy', overwrite: true

  input:
  each mutation_type from params.mutation_types
  each dataset from params.dataset
  file inputs from converted_SP_to_MSA_for_unoptimised_NNLS
  file signatures from signatures_for_unoptimised_NNLS

  output:
  file("./${dataset}/output_${dataset}_${mutation_type}_mutations_table.csv") into unoptimised_outputs
  file("./${dataset}/output_${dataset}_${mutation_type}_weights_table.csv")
  file("./${dataset}/output_${dataset}_${mutation_type}_stat_info.csv")
  set dataset, mutation_type into unoptimised_attributions

  script:
  """
  python $baseDir/bin/run_NNLS.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} -n ${params.number_of_samples} \
                                  -p ${signature_prefix} -i ${params.input_tables} -s ${params.signature_tables} -o "./"
  """
}

process run_simulations {
  tag "${mutation_type}/${dataset}"
  publishDir "$baseDir/output_tables", mode: 'copy', overwrite: true

  input:
  set dataset, mutation_type from unoptimised_attributions

  output:
  file("./SIM_${dataset}/WGS_SIM_${dataset}.${params.SBS_context}.csv") optional true
  file("./SIM_${dataset}/WGS_SIM_${dataset}.${mutation_type}.csv") optional true
  file("./SIM_${dataset}/WGS_SIM_${dataset}.dinucs.csv") optional true
  file("./SIM_${dataset}/WGS_SIM_${dataset}.indels.csv") optional true
  file("./SIM_${dataset}/WGS_SIM_${dataset}.${params.SBS_context}.weights.csv") optional true
  file("./SIM_${dataset}/WGS_SIM_${dataset}.${mutation_type}.weights.csv") optional true
  file("./SIM_${dataset}/WGS_SIM_${dataset}.dinucs.weights.csv") optional true
  file("./SIM_${dataset}/WGS_SIM_${dataset}.indels.weights.csv") optional true
  set dataset, mutation_type into simulation_outputs
  set dataset, mutation_type into simulation_outputs_for_bootstrap

  script:
  """
  python $baseDir/bin/simulate_data.py -d SIM_${dataset} -t ${mutation_type} -c ${params.SBS_context} -n ${number_of_simulated_samples} -p ${signature_prefix} \
                                  -B ${noise_flag} --noise_type ${params.noise_type} -Z ${params.noise_stdev} --zero_inflation_threshold ${params.zero_inflation_threshold} \
                                  -i $baseDir/output_tables_unoptimised/${dataset}/output_${dataset}_${mutation_type}_mutations_table.csv -s ${params.signature_tables} -o "./"
  """
}

process run_optimisation_NNLS {
  tag "${mutation_type}/${dataset}"
  publishDir "${params.optimisation_NNLS_output_path}"

  input:
  each weak_threshold from weak_thresholds
  each strong_threshold from strong_thresholds
  set dataset, mutation_type from simulation_outputs

  output:
  set dataset, mutation_type into optimisation_attribution_for_tables
  file("SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/*.csv")

  when:
  !params.run_only_simulations

  script:
  """
  mkdir -p ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}
  [[ ${mutation_type} == "SBS" ]] && \
    cp $baseDir/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.${params.SBS_context}.weights.csv ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
  [[ ${mutation_type} == "DBS" ]] && \
    cp $baseDir/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.dinucs.weights.csv ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
  [[ ${mutation_type} == "ID" ]] && \
    cp $baseDir/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.indels.weights.csv ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
  [[ ${mutation_type} == "SV" || ${mutation_type} == "CNV" ]] && \
    cp $baseDir/output_tables/SIM_${dataset}/WGS_SIM_${dataset}.${mutation_type}.weights.csv ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/
  python $baseDir/bin/run_NNLS.py -d SIM_${dataset} -t ${mutation_type} -p ${signature_prefix} \
                                  --optimisation_strategy ${params.optimisation_strategy} \
                                  -W ${weak_threshold} -S ${strong_threshold} \
                                  -i $baseDir/output_tables -s ${params.signature_tables} \
                                  -o "./" -x -c ${params.SBS_context} --add_suffix
  """
}

process run_optimisation_NNLS_bootstrapping {
  tag "${mutation_type}/${dataset}"
  publishDir "${params.optimisation_NNLS_output_path}"

  input:
  each i from 1..number_of_bootstrapped_samples_in_optimisation
  each weak_threshold from weak_thresholds
  each strong_threshold from strong_thresholds
  set dataset, mutation_type from simulation_outputs_for_bootstrap

  output:
  file("./SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output/output_SIM_${dataset}_${mutation_type}_${weak_threshold}_${strong_threshold}_${i}_mutations_table.csv")
  file("./SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output/output_SIM_${dataset}_${mutation_type}_${weak_threshold}_${strong_threshold}_${i}_stat_info.csv")
  file("./SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output/output_SIM_${dataset}_${mutation_type}_${weak_threshold}_${strong_threshold}_${i}_weights_table.csv")
  val i into optimisation_bootstrap_outputs

  when:
  !params.run_only_simulations

  script:
  """
  python $baseDir/bin/run_NNLS.py -B -d SIM_${dataset} -t ${mutation_type} -c ${params.SBS_context} -x \
                                  --optimisation_strategy ${params.optimisation_strategy} \
                                  --bootstrap_method ${params.bootstrap_method} \
                                  -W ${weak_threshold} -S ${strong_threshold} --add_suffix \
                                  -p ${signature_prefix} -i $baseDir/output_tables -s ${params.signature_tables} -o "./"
  mkdir -p SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output
  mv SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/output_SIM_${dataset}_${mutation_type}_mutations_table.csv SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output/output_SIM_${dataset}_${mutation_type}_${weak_threshold}_${strong_threshold}_${i}_mutations_table.csv
  mv SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/output_SIM_${dataset}_${mutation_type}_weights_table.csv SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output/output_SIM_${dataset}_${mutation_type}_${weak_threshold}_${strong_threshold}_${i}_weights_table.csv
  mv SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/output_SIM_${dataset}_${mutation_type}_stat_info.csv SIM_${dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${strong_threshold}/bootstrap_output/output_SIM_${dataset}_${mutation_type}_${weak_threshold}_${strong_threshold}_${i}_stat_info.csv
  """
}

process make_optimisation_bootstrap_tables {
  tag "${mutation_type}/${dataset}"
  publishDir "${params.optimisation_NNLS_output_path}"

  input:
  set dataset, mutation_type from optimisation_attribution_for_tables
  val i from optimisation_bootstrap_outputs.collect()
  each weak_threshold from weak_thresholds
  each strong_threshold from strong_thresholds

  output:
  file("*/*.csv") into all_optimisation_bootstrap_outputs_for_penalties
  file '*/truth_studies/*.csv'
  file '*/truth_studies/*.json'

  script:
  """
  python $baseDir/bin/make_bootstrap_tables.py -d SIM_${dataset} -t ${mutation_type} -p ${signature_prefix} \
          --suffix ${weak_threshold}_${strong_threshold} -l ${params.confidence_level} \
          -c ${params.SBS_context} -S ${params.signature_tables} \
          -T ${params.signature_attribution_thresholds.join(' ')} \
          -i ${params.optimisation_NNLS_output_path} -o "./" -n ${number_of_bootstrapped_samples_in_optimisation}
  """
}

process calculate_optimal_penalties {
  tag "${mutation_type}/${dataset}"
  publishDir "${params.optimisation_NNLS_output_path}"

  input:
  file ('*.csv') from all_optimisation_bootstrap_outputs_for_penalties.collect()
  each mutation_type from params.mutation_types
  each dataset from params.dataset

  output:
  set dataset, mutation_type into penalties_for_optimisation_plotting
  set dataset, mutation_type into penalties_for_central_NNLS_attribution
  set dataset, mutation_type into penalties_for_bootstrap_NNLS_attribution
  file '*/*/*.csv'
  file '*/*/*.json'
  file("SIM_${dataset}/${mutation_type}/optimal_weak_penalty") into optimal_weak_penalty
  file("SIM_${dataset}/${mutation_type}/optimal_weak_penalty") into optimal_weak_penalty_for_bootstrapping
  file("SIM_${dataset}/${mutation_type}/optimal_strong_penalty") into optimal_strong_penalty
  file("SIM_${dataset}/${mutation_type}/optimal_strong_penalty") into optimal_strong_penalty_for_bootstrapping

  script:
  """
  python $baseDir/bin/calculate_optimal_penalties.py -d SIM_${dataset} -t ${mutation_type} \
                                                     -I $baseDir/output_tables/SIM_${dataset} \
                                                     -i ${params.optimisation_NNLS_output_path} -o "./" \
                                                     -c ${params.SBS_context} \
                                                     ${no_CI_for_penalties_flag} \
                                                     ${calculate_penalty_on_average_flag} \
                                                     ${prioritised_signatures_flag} \
                                                     -M ${params.metric_to_prioritise} \
                                                     -T ${params.metric_threshold} \
                                                     -W ${weak_thresholds.join(' ')} \
                                                     -S ${strong_thresholds.join(' ')} \
                                                     --signature_path ${params.signature_tables} \
                                                     -p ${signature_prefix}
  """
}

process plot_optimisation_plots {
  tag "${mutation_type}/${dataset}"
  publishDir "${params.optimisation_plots_output_path}"

  input:
  set dataset, mutation_type from penalties_for_optimisation_plotting

  when:
  params.plot_optimisation_plots

  output:
  file '*/*.pdf' optional true
  file '*/*/*.pdf' optional true
  file '*/*/*/*.pdf' optional true

  script:
  """
  python $baseDir/bin/plot_metric_heatmaps.py -d SIM_${dataset} -t ${mutation_type} \
                                              -i ${params.optimisation_NNLS_output_path} -o "./" \
                                              -c ${params.SBS_context} \
                                              -l ${params.metric_threshold} \
                                              -W ${weak_thresholds.join(' ')} \
                                              -S ${strong_thresholds.join(' ')} \
                                              -T ${params.signature_attribution_thresholds.join(' ')} \
                                              --signature_path ${params.signature_tables} \
                                              -p ${signature_prefix}
  """
}


process run_NNLS_normal {
  tag "${mutation_type}/${dataset}"
  publishDir "$baseDir/output_tables", mode: 'copy', overwrite: true

  input:
  set dataset, mutation_type from penalties_for_central_NNLS_attribution
  file weak_penalty from optimal_weak_penalty
  file strong_penalty from optimal_strong_penalty

  output:
  file("./${dataset}/output_${dataset}_${mutation_type}_mutations_table.csv") into final_outputs_no_bootstrap
  file("./${dataset}/output_${dataset}_${mutation_type}_weights_table.csv")
  file("./${dataset}/output_${dataset}_${mutation_type}_stat_info.csv")
  file("./${dataset}/output_${dataset}_${mutation_type}_fitted_values.csv")
  file("./${dataset}/output_${dataset}_${mutation_type}_residuals.csv") into central_NNLS_residuals
  set dataset, mutation_type into attribution_for_bootstrap_plots
  set dataset, mutation_type into attribution_for_spectra_plots
  set dataset, mutation_type into attribution_for_residuals
  set dataset, mutation_type into attribution_for_metrics
  set dataset, mutation_type into attribution_for_tables

  when:
  !params.run_only_optimisation

  script:
  """
  python $baseDir/bin/run_NNLS.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} ${optimised_flag} \
                                  --optimisation_strategy ${params.optimisation_strategy} \
                                  -W `< ${weak_penalty}` -S `< ${strong_penalty}` -n ${params.number_of_samples} \
                                  -p ${signature_prefix} -i ${params.input_tables} -s ${params.signature_tables} -o "./"
  cp ${dataset}/output_${dataset}_${mutation_type}_residuals.csv ${params.input_tables}/${dataset}/
  cp ${dataset}/output_${dataset}_${mutation_type}_fitted_values.csv ${params.input_tables}/${dataset}/
  """
}

process run_NNLS_bootstrapping {
  tag "${mutation_type}/${dataset}"
  publishDir "$baseDir/output_tables", mode: 'copy', overwrite: true

  input:
  each i from 1..number_of_bootstrapped_samples
  set dataset, mutation_type from penalties_for_bootstrap_NNLS_attribution
  file weak_penalty from optimal_weak_penalty_for_bootstrapping
  file strong_penalty from optimal_strong_penalty_for_bootstrapping
  // file residuals from central_NNLS_residuals // uncomment in using residuals bootstrapping

  output:
  file("./${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_${i}_mutations_table.csv")
  file("./${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_${i}_stat_info.csv")
  file("./${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_${i}_weights_table.csv")
  val i into bootstrap_outputs

  when:
  !params.run_only_optimisation

  script:
  """
  python $baseDir/bin/run_NNLS.py -B -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} ${optimised_flag} \
                                  --optimisation_strategy ${params.optimisation_strategy} --bootstrap_method ${params.bootstrap_method} \
                                  -W `< ${weak_penalty}` -S `< ${strong_penalty}` -n ${params.number_of_samples} \
                                  -p ${signature_prefix} -i ${params.input_tables} -s ${params.signature_tables} -o "./"
  mkdir -p ${dataset}/bootstrap_output
  mv ${dataset}/output_${dataset}_${mutation_type}_mutations_table.csv ${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_${i}_mutations_table.csv
  mv ${dataset}/output_${dataset}_${mutation_type}_weights_table.csv ${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_${i}_weights_table.csv
  mv ${dataset}/output_${dataset}_${mutation_type}_stat_info.csv ${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_${i}_stat_info.csv
  """
}

process make_bootstrap_tables {
  tag "${mutation_type}/${dataset}"
  publishDir "$baseDir/output_tables", mode: 'copy', overwrite: true

  input:
  set dataset, mutation_type from attribution_for_tables
  val i from bootstrap_outputs.collect()

  output:
  file("./${dataset}/CIs_${dataset}_${mutation_type}_bootstrap_output_${suffix}.csv") into final_outputs_post_bootstrap
  file("./${dataset}/signatures_prevalences_${dataset}_${mutation_type}.csv") into signature_prevalences
  file("./${dataset}/attributions_per_sample_${dataset}_${mutation_type}_bootstrap_output_${suffix}.json") into attributions_per_sample
  file("./${dataset}/attributions_per_signature_${dataset}_${mutation_type}_bootstrap_output_${suffix}.json")
  file("./${dataset}/stat_metrics_${dataset}_${mutation_type}_bootstrap_output_${suffix}.json")
  file("./${dataset}/pruned_attribution_${dataset}_${mutation_type}_abs_mutations.csv")
  file '*/truth_studies/*.csv' optional true
  file '*/truth_studies/*.json' optional true

  script:
  """
  mkdir -p $baseDir/output_tables/${dataset}
  [[ ${dataset} == *"SIM"* ]] && [[ ${mutation_type} == "SBS" ]] && \
    cp ${params.input_tables}/${dataset}/WGS_${dataset}.${params.SBS_context}.weights.csv $baseDir/output_tables/${dataset}/
  [[ ${dataset} == *"SIM"* ]] && [[ ${mutation_type} == "DBS" ]] && \
    cp ${params.input_tables}/${dataset}/WGS_${dataset}.dinucs.weights.csv $baseDir/output_tables/${dataset}/
  [[ ${dataset} == *"SIM"* ]] && [[ ${mutation_type} == "ID" ]] && \
    cp ${params.input_tables}/${dataset}/WGS_${dataset}.indels.weights.csv $baseDir/output_tables/${dataset}/
  [[ ${dataset} == *"SIM"* ]] && [[ ${mutation_type} == "SV" || ${mutation_type} == "CNV" ]] && \
    cp ${params.input_tables}/${dataset}/WGS_${dataset}.${mutation_type}.weights.csv $baseDir/output_tables/${dataset}/
  python $baseDir/bin/make_bootstrap_tables.py -d ${dataset} -t ${mutation_type} -p ${signature_prefix} ${abs_flag} \
                                               -c ${params.SBS_context} -S ${params.signature_tables} -l ${params.confidence_level} \
                                               -T ${params.signature_attribution_thresholds.join(' ')} \
                                               -i $baseDir/output_tables -o "./" -n ${number_of_bootstrapped_samples}
  """
}

process plot_bootstrap_attributions {
  tag "${mutation_type}/${dataset}"
  publishDir "${params.plots_output_path}", mode: 'move', overwrite: true

  input:
  set dataset, mutation_type from attribution_for_bootstrap_plots
  file bootstrap_attributions from attributions_per_sample.collect()

  when:
  params.plot_bootstrap_attributions

  output:
  file '*/*/bootstrap_plots/*.pdf' optional true
  file '*/*/bootstrap_plots/*/*.pdf' optional true
  file '*/*/bootstrap_plots/*/*/*.pdf' optional true

  script:
  """
  python $baseDir/bin/plot_bootstrap_attributions.py -d ${dataset} -t ${mutation_type} -p ${signature_prefix} ${abs_flag} \
                                                     -c ${params.SBS_context} -S ${params.signature_tables} -I ${params.input_tables} \
                                                     -i $baseDir/output_tables -o "./" -n ${number_of_bootstrapped_samples}
  """
}

process plot_metrics {
  tag "${mutation_type}/${dataset}"
  publishDir "${params.plots_output_path}", mode: 'move', overwrite: true

  input:
  set dataset, mutation_type from attribution_for_metrics
  file prevalences from signature_prevalences.collect()

  when:
  params.plot_metrics

  output:
  file '*/*/bootstrap_plots/*.pdf' optional true
  file '*/*/bootstrap_plots/*/*.pdf' optional true
  file '*/*/bootstrap_plots/*/*/*.pdf' optional true

  script:
  """
  python $baseDir/bin/plot_metrics.py -d ${dataset} -t ${mutation_type} -l ${params.metric_threshold} -i $baseDir/output_tables -o "./"
  """
}

process plot_fitted_spectra {
  tag "${mutation_type}/${dataset}"
  publishDir "${params.plots_output_path}", mode: 'move'

  input:
  set dataset, mutation_type from attribution_for_spectra_plots

  output:
  file '*/*/*.pdf' optional true
  file '*/*/*/*.pdf' optional true
  file '*/*/*/*/*.pdf' optional true
  file '*/*/*/*/*/*.pdf' optional true

  when:
  params.plot_fitted_spectra

  script:
  """
  python $baseDir/bin/plot_mutation_spectra.py -f -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                -i $baseDir/output_tables ${error_flag} ${strands_flag} ${nontranscribed_flag} -o "./"
  python $baseDir/bin/plot_mutation_spectra.py -f -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                -r -i $baseDir/output_tables ${error_flag} ${strands_flag} ${nontranscribed_flag} -o "./"
  python $baseDir/bin/plot_mutation_spectra.py -C -f -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                -i $baseDir/output_tables ${error_flag} ${strands_flag} ${nontranscribed_flag} -o "./"
  python $baseDir/bin/plot_mutation_spectra.py -C -f -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                -r -i $baseDir/output_tables ${error_flag} ${strands_flag} ${nontranscribed_flag} -o "./"
  """
}

process plot_residuals {
  tag "${mutation_type}/${dataset}"
  publishDir "${params.plots_output_path}", mode: 'move'

  input:
  set dataset, mutation_type from attribution_for_residuals

  output:
  file '*/*/*.pdf' optional true
  file '*/*/*/*.pdf' optional true
  file '*/*/*/*/*.pdf' optional true
  file '*/*/*/*/*/*.pdf' optional true

  when:
  params.plot_residuals

  script:
  """
  python $baseDir/bin/plot_mutation_spectra.py -H -R -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                -i $baseDir/output_tables ${error_flag} ${strands_flag} ${nontranscribed_flag} -o "./"
  python $baseDir/bin/plot_mutation_spectra.py -C -R -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                -i $baseDir/output_tables ${error_flag} ${strands_flag} ${nontranscribed_flag} -o "./"
  """
}

process move_bootstrap_outputs {
  tag "${dataset}"
  publishDir "${params.tables_output_path}", mode: 'move', overwrite: true

  input:
  file output from final_outputs_post_bootstrap.collect()
  each dataset from params.dataset
  path output_path from "$baseDir/output_tables"

  output:
  file("MSA_output_${dataset}.tar.gz")

  shell:
  """
  cp -r ${output_path}/${dataset} .
  mkdir outputs_optimisation
  cp -r ${params.optimisation_NNLS_output_path}/SIM_${dataset} outputs_optimisation/
  cp -r ${params.optimisation_NNLS_output_path}/SIM_${dataset}_${params.SBS_context}* outputs_optimisation
  tar cvfh MSA_output_${dataset}.tar.gz ${dataset} outputs_optimisation
  rm -rf ${dataset}
  rm -rf outputs_optimisation
  """
}
