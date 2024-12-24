""" make_bootstrap_tables.py
Module make the bootstrap tables as well as calculate various metrics and scores in case of simulated data
"""

from argparse import ArgumentParser
import copy
import numpy as np
import pandas as pd
from datetime import datetime
from common_methods import make_folder_if_not_exists, write_data_to_JSON, calculate_confidence_interval, calculate_stat_scores, calculate_sensitivity_CI, calculate_specificity_CI

def calculate_sensitivity_thresholds(signatures_CPs_dict, signatures_to_consider, signature_attribution_thresholds):
    global confidence_level
    sensitivity_thresholds = pd.DataFrame(columns=signatures_to_consider, dtype=float)
    for signature in signatures_to_consider:
        lowest_threshold_scanned = signature_attribution_thresholds[0]
        highest_threshold_scanned = signature_attribution_thresholds[-1]
        if signatures_CPs_dict[lowest_threshold_scanned].loc[0, signature] >= confidence_level:
            sensitivity_thresholds.loc[0, signature] = lowest_threshold_scanned
            continue
        elif signatures_CPs_dict[highest_threshold_scanned].loc[0, signature] < confidence_level:
            sensitivity_thresholds.loc[0, signature] = np.nan
            continue
        else:
            for threshold_id, threshold in enumerate(signature_attribution_thresholds):
                if threshold_id != len(signature_attribution_thresholds)-1:
                    next_threshold = signature_attribution_thresholds[threshold_id+1]
                else:
                    next_threshold = threshold
                lower_CP = signatures_CPs_dict[threshold].loc[0, signature]
                higher_CP = signatures_CPs_dict[next_threshold].loc[0, signature]
                if lower_CP < confidence_level and higher_CP >= confidence_level:
                    sensitivity_thresholds.loc[0, signature] = threshold + (confidence_level-lower_CP)*(next_threshold-threshold)/(higher_CP-lower_CP)
                    continue
                elif lower_CP==confidence_level and higher_CP==confidence_level:
                    sensitivity_thresholds.loc[0, signature] = threshold
                    continue
    return sensitivity_thresholds

if __name__ == '__main__':
    start_time = datetime.now()
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_attributions_folder", dest="input_attributions_folder", default='output_tables/',
                        help="set path to NNLS output data")
    parser.add_argument("-S", "--signature_path", dest="signature_tables_path", default='signature_tables/',
                        help="set path to signature tables")
    parser.add_argument("-p", "--signature_prefix", dest="signatures_prefix", default='sigProfiler',
                        help="set prefix in signature filenames (sigProfiler by default)")
    parser.add_argument("-d", "--dataset", dest="dataset_name", default='SIM',
                        help="set the dataset name (e.g. SIM)")
    parser.add_argument("-o", "--output_folder", dest="output_folder", default='output_tables/',
                        help="set path to save plots")
    parser.add_argument("-t", "--mutation_type", dest="mutation_type", default='',
                        help="set mutation type (SBS, DBS, ID, SV, CNV)")
    parser.add_argument("-c", "--context", dest="context", default=192, type=int,
                        help="set SBS context (96, 192, 288)")
    parser.add_argument("-l", "--confidence_level", dest="confidence_level", default=0.95, type=float,
                        help="specify the confidence level for CI calculation (default: 0.95)")
    parser.add_argument("-a", "--use_absolute_numbers", dest="abs_numbers", action="store_true",
                        help="consider absolute numbers of mutations (relative by default)")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                        help="print additional information for debugging")
    parser.add_argument("-n", "--number_of_b_samples", dest="number_of_b_samples", default=1000, type=int,
                        help="Number of bootstrapped samples (1000 by default)")
    parser.add_argument("-N", "--number_of_samples", dest="number_of_samples", default=-1, type=int,
                        help="limit the number of samples to analyse (all by default)")
    parser.add_argument("-T", "--signature_attribution_thresholds", nargs='+', dest="signature_attribution_thresholds",
                        help="set the list of signature attribution thresholds for confidence probability plots (truth study)")
    parser.add_argument("--suffix", dest="suffix", default='',
                        help="add suffixes to inputs (useful during optimisation scan analysis)")

    args = parser.parse_args()

    dataset_name = args.dataset_name
    mutation_type = args.mutation_type
    context = args.context
    confidence_level = args.confidence_level
    number_of_b_samples = args.number_of_b_samples
    signature_tables_path = args.signature_tables_path
    signatures_prefix = args.signatures_prefix
    if not args.suffix:
        input_attributions_folder = args.input_attributions_folder + '/' + dataset_name + '/'
        output_folder = args.output_folder + '/' + dataset_name
    else:
        input_attributions_folder = args.input_attributions_folder + '/' + dataset_name + '_' + str(context) + '_NNLS_' + args.suffix + '/'
        output_folder = args.output_folder + '/' + dataset_name + '_' + str(context) + '_NNLS_' + args.suffix + '/'
    make_folder_if_not_exists(output_folder)
    if 'SIM' in dataset_name:
        make_folder_if_not_exists(output_folder + '/truth_studies')

    if not mutation_type:
        parser.error("Please specify the mutation type using -t option, e.g. add '-t SBS' to the command (DBS, ID).")
    elif mutation_type not in ['SBS', 'DBS', 'ID', 'SV', 'CNV']:
        raise ValueError("Unknown mutation type: %s. Known types: SBS, DBS, ID, SV, CNV" % mutation_type)

    print("*"*50)
    print("Making tables for %s mutation type, %s dataset" % (mutation_type, dataset_name))

    central_attribution_table_abs = pd.read_csv(input_attributions_folder + '/output_%s_%s_mutations_table.csv' % (dataset_name, mutation_type), index_col=0)
    central_attribution_table_weights = pd.read_csv(input_attributions_folder + '/output_%s_%s_weights_table.csv' % (dataset_name, mutation_type), index_col=0)
    central_stat_table = pd.read_csv(input_attributions_folder + '/output_%s_%s_stat_info.csv' % (dataset_name, mutation_type), index_col=0)
    bootstrap_attribution_table_abs_filename = input_attributions_folder + '/bootstrap_output/output_%s_%s_i_mutations_table.csv' % (dataset_name, mutation_type)
    bootstrap_attribution_table_weights_filename = input_attributions_folder + '/bootstrap_output/output_%s_%s_i_weights_table.csv' % (dataset_name, mutation_type)
    bootstrap_stat_table_filename = input_attributions_folder + '/bootstrap_output/output_%s_%s_i_stat_info.csv' % (dataset_name, mutation_type)
    if args.suffix:
        bootstrap_attribution_table_abs_filename = input_attributions_folder + '/bootstrap_output/output_%s_%s_%s_i_mutations_table.csv' % (dataset_name, mutation_type, args.suffix)
        bootstrap_attribution_table_weights_filename = input_attributions_folder + '/bootstrap_output/output_%s_%s_%s_i_weights_table.csv' % (dataset_name, mutation_type, args.suffix)
        bootstrap_stat_table_filename = input_attributions_folder + '/bootstrap_output/output_%s_%s_%s_i_stat_info.csv' % (dataset_name, mutation_type, args.suffix)

    if mutation_type == 'SBS':
        if context == 96:
            signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=[0, 1])
        elif context in [192, 288]:
            signatures = pd.read_csv('%s/%s_%s_%i_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type, context), index_col=[0, 1, 2])
        elif context == 1536:
            signatures = pd.read_csv('%s/%s_%s_%i_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type, context), index_col=0)
        else:
            raise ValueError("Context %i is not supported." % context)
    else:
        signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=0)

    if 'SIM' in dataset_name:
        if mutation_type == 'SBS':
            truth_attribution_table = pd.read_csv(input_attributions_folder + '/WGS_%s.%i.weights.csv' % (dataset_name, context), index_col=0)
        elif mutation_type == 'DBS':
            truth_attribution_table = pd.read_csv(input_attributions_folder + '/WGS_%s.dinucs.weights.csv' % dataset_name, index_col=0)
        elif mutation_type == 'ID':
            truth_attribution_table = pd.read_csv(input_attributions_folder + '/WGS_%s.indels.weights.csv' % dataset_name, index_col=0)
        else: # SV and CNV
            truth_attribution_table = pd.read_csv(input_attributions_folder + '/WGS_%s.%s.weights.csv' % (dataset_name, mutation_type), index_col=0)

    if args.abs_numbers:
        central_attribution_table = central_attribution_table_abs
        bootstrap_attribution_table_filename = bootstrap_attribution_table_abs_filename
        filename = dataset_name + '_' + mutation_type + '_bootstrap_output_abs_mutations'
    else:
        central_attribution_table = central_attribution_table_weights
        bootstrap_attribution_table_filename = bootstrap_attribution_table_weights_filename
        filename = dataset_name + '_' + mutation_type + '_bootstrap_output_weights'

    # limit the number of samples to analyse (if specified by -N option)
    if args.number_of_samples != -1:
        central_attribution_table = central_attribution_table.head(args.number_of_samples)
        if 'SIM' in dataset_name:
            truth_attribution_table = truth_attribution_table.head(args.number_of_samples)

    # # convert columns and index to str
    # central_attribution_table.columns = central_attribution_table.columns.astype(str)
    # central_attribution_table.index = central_attribution_table.index.astype(str)
    # if 'SIM' in dataset_name:
    #     truth_attribution_table.columns = truth_attribution_table.columns.astype(str)
    #     truth_attribution_table.index = truth_attribution_table.index.astype(str)

    main_title = dataset_name.replace('_', '/') + ' data, ' + mutation_type + ' mutation type'

    # samples, signatures and metrics to consider
    samples = central_attribution_table.index.to_list()
    signatures_to_consider = list(central_attribution_table_abs.columns)
    signatures = signatures[signatures_to_consider]
    stat_metrics = central_stat_table.columns.to_list()

    # initialise statistical metrics to fill with bootstrap values
    stat_metrics_dict = {}
    for metric in stat_metrics:
        stat_metrics_dict[metric] = pd.DataFrame(index=range(number_of_b_samples), columns=samples, dtype=float)

    # initialise per-sample attribution dictionary
    attributions_per_sample_dict = {}
    for sample in samples:
        attributions_per_sample_dict[sample] = pd.DataFrame(index=range(number_of_b_samples), columns=signatures_to_consider, dtype=float)

    # initialise per-signature attribution dictionary
    attributions_per_signature_dict = {}
    for signature in signatures_to_consider:
        attributions_per_signature_dict[signature] = pd.DataFrame(index=range(number_of_b_samples), columns=samples, dtype=float)

    # mutation categories from signatures table
    categories = signatures.index.to_list()
    # initialise mutation spectra dictionary
    # mutation_spectra_dict = {}
    # for sample in samples:
    #     mutation_spectra_dict[sample] = pd.DataFrame(index=range(number_of_b_samples), columns=categories, dtype=float)

    # if SIM (synthetic) dataset, determine acting signatures
    if 'SIM' in dataset_name:
        acting_signatures = []
        for signature in list(truth_attribution_table.columns):
            if truth_attribution_table[signature].max() > 0:
                acting_signatures.append(signature)

    # initialise more dataframes and dictionaries
    signatures_prevalences = pd.DataFrame(index=central_attribution_table.index, columns=signatures_to_consider, dtype=float)
    if 'SIM' in dataset_name:
        scores = ['Sensitivity', 'Specificity', 'Precision', 'Accuracy', 'F1', 'MCC']
        stat_scores_per_sig = {}
        stat_scores_from_CI_per_sig = {}
        sensitivity_CI_per_sig = {}
        sensitivity_CI_from_CI_per_sig = {}
        specificity_CI_per_sig = {}
        specificity_CI_from_CI_per_sig = {}
        for signature in signatures_to_consider:
            stat_scores_per_sig[signature] = pd.DataFrame(columns=scores, dtype=float)
            stat_scores_from_CI_per_sig[signature] = pd.DataFrame(columns=scores, dtype=float)
            sensitivity_CI_per_sig[signature] = pd.DataFrame(columns=['lower_CL','upper_CL'], dtype=float)
            sensitivity_CI_from_CI_per_sig[signature] = pd.DataFrame(columns=['lower_CL','upper_CL'], dtype=float)
            specificity_CI_per_sig[signature] = pd.DataFrame(columns=['lower_CL','upper_CL'], dtype=float)
            specificity_CI_from_CI_per_sig[signature] = pd.DataFrame(columns=['lower_CL','upper_CL'], dtype=float)
        stat_scores_tables = pd.DataFrame(columns=scores, dtype=float)
        stat_scores_from_CI_tables = pd.DataFrame(columns=scores, dtype=float)
        sensitivity_CI_tables = pd.DataFrame(columns=['lower_CL','upper_CL'], dtype=float)
        specificity_CI_tables = pd.DataFrame(columns=['lower_CL','upper_CL'], dtype=float)
        sensitivity_CI_from_CI_tables = pd.DataFrame(columns=['lower_CL','upper_CL'], dtype=float)
        specificity_CI_from_CI_tables = pd.DataFrame(columns=['lower_CL','upper_CL'], dtype=float)
        lower_CI_attributions = pd.DataFrame(index=central_attribution_table.index, columns=signatures_to_consider, dtype=float)
        signatures_CPs = pd.DataFrame(index=central_attribution_table.index, columns=acting_signatures, dtype=float)
        truth_and_measured_difference = pd.DataFrame(0, index=central_attribution_table.index, columns=['Truth - mean', 'Truth - median', 'Truth - NNLS'], dtype=float)
        signatures_scores = {}
        signatures_scores_from_CI = {}
        for score in scores:
            signatures_scores[score] = pd.DataFrame(columns=acting_signatures, dtype=float)
            signatures_scores_from_CI[score] = pd.DataFrame(columns=acting_signatures, dtype=float)

    print('Initialisation done.')
    print('Elapsed time:', datetime.now() - start_time)

    # fill bootstrap dataframes
    for i in range(number_of_b_samples):
        bootstrap_attribution_table = pd.read_csv(bootstrap_attribution_table_filename.replace('_i_', '_%i_' % (i+1)), index_col=0)
        bootstrap_attribution_table_abs = pd.read_csv(bootstrap_attribution_table_abs_filename.replace('_i_', '_%i_' % (i+1)), index_col=0)
        stat_table = pd.read_csv(bootstrap_stat_table_filename.replace('_i_', '_%i_' % (i+1)), index_col=0)
        # back-calculate mutation spectra from bootstrap attributions and signatures
        # mutation_spectra = bootstrap_attribution_table_abs.dot(signatures.T)
        for sample in samples:
            for signature in signatures_to_consider:
                attributions_per_sample_dict[sample].loc[i, signature] = bootstrap_attribution_table.loc[sample, signature]
                attributions_per_signature_dict[signature].loc[i, sample] = bootstrap_attribution_table.loc[sample, signature]
            # for category in categories:
            #     mutation_spectra_dict[sample].loc[i, category] = mutation_spectra.loc[sample, category]
            for metric in stat_metrics:
                stat_metrics_dict[metric].loc[i, sample] = stat_table.loc[sample, metric]
                stat_metrics_dict[metric].loc[i, sample] = stat_table.loc[sample, metric]

    print('Bootstrap dataframes filled.')
    print('Elapsed time:', datetime.now() - start_time)

    # calculate confidence intervals and 'pruned' attributions, with 0 values for signatures where CIs consistent with 0
    confidence_intervals = pd.DataFrame(index=samples, columns=signatures_to_consider, dtype=object)
    pruned_attribution_table = pd.DataFrame(index=samples, columns=signatures_to_consider, dtype=int)
    for sample in samples:
        for signature in signatures_to_consider:
            confidence_interval = calculate_confidence_interval(attributions_per_sample_dict[sample][signature], confidence = confidence_level*100)
            central_value = central_attribution_table.loc[sample, signature]
            confidence_intervals.loc[sample, signature] = [central_value, confidence_interval]
            if confidence_interval[0]==0:
                pruned_attribution_table.loc[sample, signature] = 0
            else:
                pruned_attribution_table.loc[sample, signature] = np.around(central_attribution_table_abs.loc[sample, signature])

    print('Confidence intervals calculated.')
    print('Elapsed time:', datetime.now() - start_time)

    # calculate signature confidence probabilities
    if 'SIM' in dataset_name:
        for sample in samples:
            for signature in signatures_to_consider:
                # if signature is not modelled, add it as zeros in truth table:
                if signature not in truth_attribution_table.columns:
                    truth_attribution_table[signature] = 0
                array = attributions_per_sample_dict[sample][signature]
                central = confidence_intervals.loc[sample, signature][0]
                lower, upper = confidence_intervals.loc[sample, signature][1]
                lower_CI_attributions.loc[sample, signature] = lower
                if signature in acting_signatures:
                    if lower <= truth_attribution_table.loc[sample, signature] <= upper:
                        signatures_CPs.loc[sample, signature] = 1
                    else:
                        signatures_CPs.loc[sample, signature] = 0
                else:
                    if lower <= truth_attribution_table.loc[sample, signature] <= upper:
                        signatures_CPs.loc[sample, 'Others'] = 1
                    else:
                        signatures_CPs.loc[sample, 'Others'] = 0

                # fill truth and measured difference dataframe
                truth_and_measured_difference.loc[sample, 'Truth - mean'] += truth_attribution_table.loc[sample, signature] - np.mean(array)
                truth_and_measured_difference.loc[sample, 'Truth - median'] += truth_attribution_table.loc[sample, signature] - np.median(array)
                truth_and_measured_difference.loc[sample, 'Truth - NNLS'] += truth_attribution_table.loc[sample, signature] - central_attribution_table.loc[sample, signature]

        # calculate out-of-the box confidence probabilities for all acting signatures
        signatures_CPs = signatures_CPs.sum(axis=0).to_frame().T.div(len(samples))
        print('Signature confidence probabilities calculated.')
        print('Elapsed time:', datetime.now() - start_time)

    # calculate CPs for a range of sensitivity thresholds for each sig
    if 'SIM' in dataset_name:
        signature_attribution_thresholds = [float(i)/100 for i in args.signature_attribution_thresholds]
        signatures_CPs_dict = {}
        for threshold in signature_attribution_thresholds:
            signatures_CPs_dict[threshold] = pd.DataFrame(index=central_attribution_table.index, columns=signatures_to_consider, dtype=float)
            for sample in samples:
                for signature in signatures_to_consider:
                    central = confidence_intervals.loc[sample, signature][0]
                    lower, upper = confidence_intervals.loc[sample, signature][1]
                    # increase upper bound of CI if less than threshold
                    if lower == 0 and upper < threshold:
                        upper = threshold
                    if lower <= truth_attribution_table.loc[sample, signature] <= upper:
                        signatures_CPs_dict[threshold].loc[sample, signature] = 1
                    elif lower > 0 and upper > 0 and truth_attribution_table.loc[sample, signature] > 0:
                        # this case includes cases when CI is not consistent with 0, signature is simulated but CI does not contain the true value
                        signatures_CPs_dict[threshold].loc[sample, signature] = 1
                    else:
                        signatures_CPs_dict[threshold].loc[sample, signature] = 0
            # normalise to calculate CPs
            signatures_CPs_dict[threshold] = signatures_CPs_dict[threshold].sum(axis=0).to_frame().T.div(len(samples))

        # extract attribution  thresholds at which CP reaches specified confidence level for each signature
        sensitivity_thresholds = calculate_sensitivity_thresholds(signatures_CPs_dict, signatures_to_consider, signature_attribution_thresholds)
        print('Signature sensitivity thresholds calculated.')
        print('Elapsed time:', datetime.now() - start_time)

    # # correct attributions per sample dict to disregard medians consistent with 0
    # for sample in samples:
    #     for signature in signatures_to_consider:
    #         array = attributions_per_sample_dict[sample][signature]
    #         # disregard attributions consistent with 0
    #         # if confidence_interval[0]==0:
    #         # disregard attribution with 0 CI median
    #         if np.median(array) == 0:
    #             attributions_per_sample_dict[sample][signature] = 0

    # calculate signature prevalences
    for sample in samples:
        signatures_prevalences.loc[sample, :] = attributions_per_sample_dict[sample].astype(bool).sum(axis=0)
        # if 'SIM' in dataset_name:
            # # alternative procedure (considering bootstrap but not CIs)
            # truth_one_sample_table = truth_attribution_table.loc[sample, :].to_frame().T
            # truth_table = truth_one_sample_table.loc[truth_one_sample_table.index.repeat(number_of_b_samples)]
            # stat_scores_tables.loc[sample, :] = calculate_stat_scores(signatures, attributions_per_sample_dict[sample], truth_table)
            # for signature in signatures_to_consider:
            #     stat_scores_per_sig[signature].loc[sample, :] = calculate_stat_scores([signature], attributions_per_sample_dict[sample], truth_table)

    signatures_prevalences = signatures_prevalences/number_of_b_samples

    print('Signature prevalences calculated.')
    print('Elapsed time:', datetime.now() - start_time)

    # calculate per-signature statistical scores
    if 'SIM' in dataset_name:
        stat_scores_tables.loc[0, :] = calculate_stat_scores(signatures, central_attribution_table, truth_attribution_table)
        stat_scores_from_CI_tables.loc[0, :] = calculate_stat_scores(signatures, lower_CI_attributions, truth_attribution_table)
        sensitivity_CI_tables.loc[0, :] = calculate_sensitivity_CI(signatures, central_attribution_table, truth_attribution_table)
        sensitivity_CI_from_CI_tables.loc[0, :] = calculate_sensitivity_CI(signatures, lower_CI_attributions, truth_attribution_table)
        specificity_CI_tables.loc[0, :] = calculate_specificity_CI(signatures, central_attribution_table, truth_attribution_table)
        specificity_CI_from_CI_tables.loc[0, :] = calculate_specificity_CI(signatures, lower_CI_attributions, truth_attribution_table)
        for signature in signatures_to_consider:
            stat_scores_per_sig[signature].loc[0, :] = calculate_stat_scores([signature], central_attribution_table, truth_attribution_table)
            stat_scores_from_CI_per_sig[signature].loc[0, :] = calculate_stat_scores([signature], lower_CI_attributions, truth_attribution_table)
            sensitivity_CI_per_sig[signature].loc[0, :] = calculate_sensitivity_CI([signature], central_attribution_table, truth_attribution_table)
            sensitivity_CI_from_CI_per_sig[signature].loc[0, :] = calculate_sensitivity_CI([signature], lower_CI_attributions, truth_attribution_table)
            specificity_CI_per_sig[signature].loc[0, :] = calculate_specificity_CI([signature], central_attribution_table, truth_attribution_table)
            specificity_CI_from_CI_per_sig[signature].loc[0, :] = calculate_specificity_CI([signature], lower_CI_attributions, truth_attribution_table)
            for score in scores:
                signatures_scores[score].loc[0, signature] = stat_scores_per_sig[signature].loc[0, score]
                signatures_scores_from_CI[score].loc[0, signature] = stat_scores_from_CI_per_sig[signature].loc[0, score]
        print('Signature stat scores calculated.')
        print('Elapsed time:', datetime.now() - start_time)

    # write the outputs to files
    if 'SIM' in dataset_name:
        truth_and_measured_difference.to_csv(output_folder + '/truth_studies/truth_and_measured_difference_' + mutation_type + '.csv')
        signatures_CPs.to_csv(output_folder + '/truth_studies/signatures_CPs_' + mutation_type + '.csv')
        sensitivity_thresholds.to_csv(output_folder + '/truth_studies/sensitivity_thresholds_' + mutation_type + '.csv')
        stat_scores_from_CI_tables.to_csv(output_folder + '/truth_studies/stat_scores_from_CI_tables_' + mutation_type + '.csv')
        stat_scores_tables.to_csv(output_folder + '/truth_studies/stat_scores_tables_' + mutation_type + '.csv')
        sensitivity_CI_from_CI_tables.to_csv(output_folder + '/truth_studies/sensitivity_CI_from_CI_tables_' + mutation_type + '.csv')
        sensitivity_CI_tables.to_csv(output_folder + '/truth_studies/sensitivity_CI_tables_' + mutation_type + '.csv')
        specificity_CI_from_CI_tables.to_csv(output_folder + '/truth_studies/specificity_CI_from_CI_tables_' + mutation_type + '.csv')
        specificity_CI_tables.to_csv(output_folder + '/truth_studies/specificity_CI_tables_' + mutation_type + '.csv')
        write_data_to_JSON(signatures_CPs_dict, output_folder + '/truth_studies/signatures_CPs_dict_' + mutation_type + '.json')
        write_data_to_JSON(signatures_scores_from_CI, output_folder + '/truth_studies/signatures_scores_from_CI_' + mutation_type + '.json')
        write_data_to_JSON(signatures_scores, output_folder + '/truth_studies/signatures_scores_' + mutation_type + '.json')
        write_data_to_JSON(stat_scores_from_CI_per_sig, output_folder + '/truth_studies/stat_scores_from_CI_per_sig_' + mutation_type + '.json')
        write_data_to_JSON(stat_scores_per_sig, output_folder + '/truth_studies/stat_scores_per_sig_' + mutation_type + '.json')
        write_data_to_JSON(sensitivity_CI_from_CI_per_sig, output_folder + '/truth_studies/sensitivity_CI_from_CI_per_sig_' + mutation_type + '.json')
        write_data_to_JSON(sensitivity_CI_per_sig, output_folder + '/truth_studies/sensitivity_CI_per_sig_' + mutation_type + '.json')
        write_data_to_JSON(specificity_CI_from_CI_per_sig, output_folder + '/truth_studies/specificity_CI_from_CI_per_sig_' + mutation_type + '.json')
        write_data_to_JSON(specificity_CI_per_sig, output_folder + '/truth_studies/specificity_CI_per_sig_' + mutation_type + '.json')
        print('Truth studies outputs written.')
        print('Elapsed time:', datetime.now() - start_time)

    # common outputs
    confidence_intervals.to_csv(output_folder + '/CIs_' + filename + '.csv')
    signatures_prevalences.to_csv(output_folder + '/signatures_prevalences_' + dataset_name + '_' + mutation_type + '.csv')
    write_data_to_JSON(attributions_per_sample_dict, output_folder + '/attributions_per_sample_' + filename + '.json')
    write_data_to_JSON(attributions_per_signature_dict, output_folder + '/attributions_per_signature_' + filename + '.json')
    write_data_to_JSON(stat_metrics_dict, output_folder + '/stat_metrics_' + filename + '.json')
    # # pd.to_json function does not fully support MultiIndex (present in signature tables), hence commented out for now
    # write_data_to_JSON(mutation_spectra_dict, output_folder + '/mutation_spectra_' + filename + '.json')
    pruned_attribution_table.to_csv(output_folder + '/pruned_attribution_' + dataset_name + '_' + mutation_type + '_abs_mutations.csv')

    print('All remaining outputs written. Done!')
    print('Elapsed time:', datetime.now() - start_time)
