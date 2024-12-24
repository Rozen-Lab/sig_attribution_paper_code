import argparse
import os
import warnings
import copy
import pandas as pd
import numpy as np
from common_methods import make_folder_if_not_exists, calculate_similarity, read_data_from_JSON, write_data_to_JSON

# non-CI approach to similarity metrics
def measure_similarity_metrics(method):
    global input_truth_path, input_reco_path
    global dataset, context, metrics, weak_thresholds, strong_thresholds, mutation_type, truth_table_filename
    # initialise dicts and arrays
    similarity_tables = {}
    similarity_uncertainty_tables = {}

    for metric in metrics:
        similarity_tables[metric] = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
        similarity_uncertainty_tables[metric] = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)

    rss_table = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
    chi2_table = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)

    # truth_table_filename = input_truth_path + '/WGS_' + dataset + '.' + str(context) + '.weights.csv'
    truth_table = pd.read_csv(truth_table_filename, index_col=0)

    for weak_threshold in weak_thresholds:
        for strong_threshold in strong_thresholds:
            reco_table_filename = input_reco_path + '/' + dataset + '_' + \
                str(context) + '_' + method + '_' + weak_threshold + \
                '_' + strong_threshold + '/output_%s_%s_weights_table.csv' % (dataset, mutation_type)
            # reco_table_filename = input_reco_path + '/output_tables_' + weak_threshold + '_' + strong_threshold + '/' + dataset + '/output_%s_SBS_weights_table.csv' % dataset
            reco_table = pd.read_csv(reco_table_filename, index_col=0)

            # remove missing sigs (some methods pre-filter input signatures)
            truth_sigs = truth_table.columns.to_list()
            reco_sigs = reco_table.columns.to_list()
            missing_sigs = set(truth_sigs) - set(reco_sigs)
            truth_table = truth_table.drop(missing_sigs, axis=1)

            # add missing true sigs as zeros (if truth table only contains non-zero sigs)
            zero_sigs = set(reco_sigs) - set(truth_sigs)
            for zero_sig in zero_sigs:
                truth_table[zero_sig] = 0

            all_signatures = list(truth_table.columns)

            # reindex so that columns ordering is the same for both reco and truth
            truth_table = truth_table.reindex(reco_table.columns, axis=1)

            # print(reco_table.loc[:, (reco_table != 0).any(axis=0)])
            # print(truth_table.loc[:, (truth_table != 0).any(axis=0)])
            if set(reco_table.columns.to_list()) != set(truth_table.columns.to_list()):
                # print(reco_table.columns.to_list())
                # print(truth_table.columns.to_list())
                raise ValueError("Incompatible truth and reco tables: check signature columns")
            if set(reco_table.index) != set(truth_table.index):
                # print(reco_table.index.to_list())
                # print(truth_table.index.to_list())
                raise ValueError("Incompatible truth and reco tables: check indexing (samples)")

            number_of_samples = len(reco_table)

            # RSS and Chi2 for NNLS from separate stat info table
            if 'NNLS' in method:
                stat_table_filename = input_reco_path + '/' + dataset + '_' + \
                    str(context) + '_' + method + '_' + weak_threshold + \
                    '_' + strong_threshold + '/output_%s_%s_stat_info.csv' % (dataset, mutation_type)
                # stat_table_filename = input_reco_path + '/output_tables_' + weak_threshold + '_' + strong_threshold + '/' + dataset + '/output_%s_SBS_stat_info.csv' % dataset
                stat_table = pd.read_csv(stat_table_filename, index_col=0)

                rss_table.loc[weak_threshold, strong_threshold] = np.mean(
                    stat_table['RSS'].head(number_of_samples).values)
                chi2_table.loc[weak_threshold, strong_threshold] = np.mean(
                    stat_table['Chi2'].head(number_of_samples).values)

            # calculate similarities for input metrics
            for metric in metrics:
                similarities = []
                for i in range(number_of_samples):
                    reco_sample = reco_table.iloc[i].values
                    truth_sample = truth_table.iloc[i].values
                    if not truth_sample.any():
                        # similarity (e.g. cosine or JS) is undefined for zero samples, so skipping
                        continue
                    similarity = calculate_similarity(reco_sample, truth_sample, metric)
                    # if np.isnan(similarity):
                    #     warnings.warn('Warning: Nan similarity in metric %s, thresholds %s/%s, sample: %s, %s, %s' % (metric, weak_threshold, strong_threshold, i, reco_sample, truth_sample))
                    similarities.append(similarity)

                similarity_tables[metric].loc[weak_threshold, strong_threshold] = np.mean(similarities)
                similarity_uncertainty_tables[metric].loc[weak_threshold, strong_threshold] = np.std(
                    similarities)
    
    return similarity_tables, similarity_uncertainty_tables, rss_table, chi2_table

def calculate_optimal_penalty(sensitivity_table, specificity_table, label = '', metric_to_prioritise='specificity', threshold=0.95):
    assert sensitivity_table.index.equals(specificity_table.index)
    assert sensitivity_table.columns.equals(specificity_table.columns)
    if np.isnan(sensitivity_table.values).all():
        warnings.warn('Warning: %s sensitivity table is full of NaN values, returning NaN optimal penalty' % label)
        return np.nan, np.nan
    if np.isnan(specificity_table.values).all():
        warnings.warn('Warning: %s specificity table is full of NaN values, returning NaN optimal penalty' % label)
        return np.nan, np.nan
    print('Calculating optimal penalties for %s, prioritising %s with %.2f threshold' % (label, metric_to_prioritise, threshold))
    optimal_penalties = []
    if metric_to_prioritise == 'specificity':
        if specificity_table.max().max() < threshold:
            warnings.warn('Warning: %s specificity never reaches %.2f, returning best alternative which is %.2f' % (label, threshold, specificity_table.max().max()))
            # returning the location of the maximum sensitivity element, out of elements with maximum specificity
            optimal_penalties = sensitivity_table[specificity_table==specificity_table.stack().max()].stack().idxmax()
        else:
            # returning the location of the maximum sensitivity element, out of elements where specificity >= threshold
            optimal_penalties = sensitivity_table[specificity_table>=threshold].stack().idxmax()
    elif metric_to_prioritise == 'sensitivity':
        if sensitivity_table.max().max() < threshold:
            warnings.warn('Warning: %s sensitivity never reaches %.2f, returning best alternative which is %.2f' % (label, threshold, sensitivity_table.max().max()))
            # returning the location of the maximum specificity element, out of elements with maximum sensitivity
            optimal_penalties = specificity_table[sensitivity_table==sensitivity_table.stack().max()].stack().idxmax()
        else:
            # returning the location of the maximum specificity element, out of elements where sensitivity >= threshold
            optimal_penalties = specificity_table[sensitivity_table>=threshold].stack().idxmax()
    else:
        raise ValueError("Unsupported metric to prioritise: %s" % metric_to_prioritise)
    if np.nan in optimal_penalties:
        warnings.warn('Warning: nan penalties in', optimal_penalties)
    print('Optimal penalties returned:', optimal_penalties)
    return optimal_penalties

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--mutation_type", dest="mutation_type", default='',
                      help="set mutation type (SBS, DBS, ID, SV, CNV)")
    parser.add_argument("-d", "--dataset", dest="dataset_name", default='SIM',
                      help="set the dataset name ('SIM' by default)")
    parser.add_argument("-c", "--context", dest="context", default=96, type=int,
                      help="set SBS context (96, 192, 288, 1536)")
    parser.add_argument("-m", "--method", dest="method", default='NNLS',
                      help="set the method name (e.g. 'NNLS')")
    parser.add_argument("-i", "--input_reco_path", dest="input_reco_path", default='output_opt_check/',
                      help="set path to input reconstructed weights tables")
    parser.add_argument("-I", "--input_truth_path", dest="input_truth_path", default='input_mutation_tables/SIM/',
                      help="set path to input truth (simulated) weights tables")
    parser.add_argument("-o", "--output_path", dest="output_path", default='output_optimisation/',
                      help="set path to save output files")
    parser.add_argument("--no_CI", dest="no_CI", action="store_true",
                      help="do not use confidence intervals to pick optimal penalties")
    parser.add_argument("-M", "--metric_to_prioritise", dest="metric_to_prioritise", default='specificity',
                      help="set a metric to prioritise (default: specificity), requiring at least the specified threshold or closest alternative")
    parser.add_argument("-T", "--metric_threshold", dest="metric_threshold", default=0.95, type=float,
                      help="specify the minimum threshold of the prioritised metric")
    parser.add_argument("--signatures_to_prioritise", nargs='+', dest="signatures_to_prioritise",
                      help="set a list of signatures to prioritise (all by default)")
    parser.add_argument("--average", dest="average", action="store_true",
                      help="apply criteria based on signatures overall (on average, less conservative), rather than maximising prioritised metric for every signature (more conservative)")
    parser.add_argument("-W", "--weak_thresholds", nargs='+', dest="weak_thresholds",
                      help="set the list of weak thresholds")
    parser.add_argument("-S", "--strong_thresholds", nargs='+', dest="strong_thresholds",
                      help="set the list of strong thresholds")
    parser.add_argument("--signature_path", dest="signature_tables_path", default='signature_tables/',
                      help="set path to signature tables")
    parser.add_argument("-p", "--signature_prefix", dest="signatures_prefix", default='sigProfiler',
                      help="set prefix in signature filenames (sigProfiler by default)")
    args = parser.parse_args()

    mutation_type = args.mutation_type
    dataset = args.dataset_name
    input_reco_path = args.input_reco_path
    input_truth_path = args.input_truth_path
    method = args.method
    context = args.context
    signature_tables_path = args.signature_tables_path
    signatures_prefix = args.signatures_prefix
    weak_thresholds = args.weak_thresholds
    strong_thresholds = args.strong_thresholds

    no_CI = args.no_CI
    metric_to_prioritise = args.metric_to_prioritise
    metric_threshold = args.metric_threshold
    average = args.average
    if args.signatures_to_prioritise:
        signatures_to_prioritise = [sig for sig in args.signatures_to_prioritise if mutation_type in sig] # only keeping relevant mutation type signatures
    else:
        signatures_to_prioritise = []

    if not weak_thresholds or not strong_thresholds:
        raise ValueError("Please specify the lists of thresholds.")

    if not mutation_type:
        parser.error("Please specify the mutation type using -t option, e.g. add '-t SBS' to the command (DBS, ID).")
    elif mutation_type not in ['SBS','DBS','ID','SV','CNV']:
        raise ValueError("Unknown mutation type: %s. Known types: SBS, DBS, ID" % mutation_type)

    metrics = ['Cosine', 'Correlation', 'L1', 'L2', 'Chebyshev', 'Jensen-Shannon']

    output_path = args.output_path + '/' + dataset + '/' + mutation_type

    make_folder_if_not_exists(output_path)

    if not method:
        parser.error("Please specify the method name (e.g. 'NNLS', 'SigProfiler') using -m option.")
    if not input_reco_path:
        parser.error("Please specify the input path for reconstructed weights tables using -i option.")
    if not input_truth_path:
        parser.error("Please specify the input truth (simulated) weights table using -I option.")

    if mutation_type=='SBS':
        truth_table_filename = input_truth_path + '/WGS_%s.%i.weights.csv' % (dataset, context)
        if context == 96:
            signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=[0, 1])
        elif context in [192, 288]:
            signatures = pd.read_csv('%s/%s_%s_%i_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type, context), index_col=[0, 1, 2])
        elif context == 1536:
            signatures = pd.read_csv('%s/%s_%s_%i_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type, context), index_col=0)
        else:
            raise ValueError("Context %i is not supported." % context)
    elif mutation_type=='DBS':
        truth_table_filename = input_truth_path + '/WGS_%s.dinucs.weights.csv' % dataset
        signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=0)
    elif mutation_type=='ID':
        truth_table_filename = input_truth_path + '/WGS_%s.indels.weights.csv' % dataset
        signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=0)
    else: # SV and CNV
        truth_table_filename = input_truth_path + '/WGS_%s.%s.weights.csv' % (dataset, mutation_type)
        signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=0)

    # truth_table = pd.read_csv(truth_table_filename, index_col=0)
    all_signatures = list(signatures.columns)
    print('Signatures analysed:', all_signatures)

    sensitivity_tables_per_sig_from_CI = {}
    specificity_tables_per_sig_from_CI = {}
    sensitivity_CI_tables_per_sig_from_CI = {}
    specificity_CI_tables_per_sig_from_CI = {}
    precision_tables_per_sig_from_CI = {}
    accuracy_tables_per_sig_from_CI = {}
    MCC_tables_per_sig_from_CI = {}
    sensitivity_tables_per_sig = {}
    specificity_tables_per_sig = {}
    sensitivity_CI_tables_per_sig = {}
    specificity_CI_tables_per_sig = {}
    precision_tables_per_sig = {}
    accuracy_tables_per_sig = {}
    MCC_tables_per_sig = {}
    sensitivity_table = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
    specificity_table = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
    sensitivity_CI_table = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
    specificity_CI_table = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
    precision_table = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
    accuracy_table = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
    MCC_table = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
    sensitivity_table_from_CI = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
    specificity_table_from_CI = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
    sensitivity_CI_table_from_CI = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
    specificity_CI_table_from_CI = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
    precision_table_from_CI = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
    accuracy_table_from_CI = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
    MCC_table_from_CI = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
    for signature in all_signatures:
        sensitivity_tables_per_sig_from_CI[signature] = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
        specificity_tables_per_sig_from_CI[signature] = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
        sensitivity_CI_tables_per_sig_from_CI[signature] = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
        specificity_CI_tables_per_sig_from_CI[signature] = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
        precision_tables_per_sig_from_CI[signature] = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
        accuracy_tables_per_sig_from_CI[signature] = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
        MCC_tables_per_sig_from_CI[signature] = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
        sensitivity_tables_per_sig[signature] = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
        specificity_tables_per_sig[signature] = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
        sensitivity_CI_tables_per_sig[signature] = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
        specificity_CI_tables_per_sig[signature] = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
        precision_tables_per_sig[signature] = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
        accuracy_tables_per_sig[signature] = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
        MCC_tables_per_sig[signature] = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)

    for weak_threshold in weak_thresholds:
        for strong_threshold in strong_thresholds:
            input_attributions_folder = input_reco_path + '/' + dataset + '_' + str(context) + '_' + method + '_' + weak_threshold + '_' + strong_threshold

            # sensitivity_thresholds = pd.read_csv(input_attributions_folder + '/truth_studies/sensitivity_thresholds_' + mutation_type + '.csv', index_col=0)
            # signatures_scores = read_data_from_JSON(input_attributions_folder + '/truth_studies/signatures_scores_' + mutation_type + '.json')
            stat_scores_tables = pd.read_csv(input_attributions_folder + '/truth_studies/stat_scores_tables_' + mutation_type + '.csv', index_col=0)
            stat_scores_from_CI_tables = pd.read_csv(input_attributions_folder + '/truth_studies/stat_scores_from_CI_tables_' + mutation_type + '.csv', index_col=0)
            stat_scores_from_CI_per_sig = read_data_from_JSON(input_attributions_folder + '/truth_studies/stat_scores_from_CI_per_sig_' + mutation_type + '.json')
            stat_scores_per_sig = read_data_from_JSON(input_attributions_folder + '/truth_studies/stat_scores_per_sig_' + mutation_type + '.json')
            sensitivity_CI_tables = pd.read_csv(input_attributions_folder + '/truth_studies/sensitivity_CI_tables_' + mutation_type + '.csv', index_col=0)
            specificity_CI_tables = pd.read_csv(input_attributions_folder + '/truth_studies/specificity_CI_tables_' + mutation_type + '.csv', index_col=0)
            sensitivity_CI_per_sig = read_data_from_JSON(input_attributions_folder + '/truth_studies/sensitivity_CI_per_sig_' + mutation_type + '.json')
            specificity_CI_per_sig = read_data_from_JSON(input_attributions_folder + '/truth_studies/specificity_CI_per_sig_' + mutation_type + '.json')
            sensitivity_CI_from_CI_tables = pd.read_csv(input_attributions_folder + '/truth_studies/sensitivity_CI_from_CI_tables_' + mutation_type + '.csv', index_col=0)
            specificity_CI_from_CI_tables = pd.read_csv(input_attributions_folder + '/truth_studies/specificity_CI_from_CI_tables_' + mutation_type + '.csv', index_col=0)
            sensitivity_CI_from_CI_per_sig = read_data_from_JSON(input_attributions_folder + '/truth_studies/sensitivity_CI_from_CI_per_sig_' + mutation_type + '.json')
            specificity_CI_from_CI_per_sig = read_data_from_JSON(input_attributions_folder + '/truth_studies/specificity_CI_from_CI_per_sig_' + mutation_type + '.json')

            sensitivity_table_from_CI.at[weak_threshold, strong_threshold] = stat_scores_from_CI_tables['Sensitivity'][0]
            specificity_table_from_CI.at[weak_threshold, strong_threshold] = stat_scores_from_CI_tables['Specificity'][0]
            sensitivity_CI_table_from_CI.at[weak_threshold, strong_threshold] = '[' + str(sensitivity_CI_from_CI_tables['lower_CL'][0]) + ', ' + str(sensitivity_CI_from_CI_tables['upper_CL'][0]) + ']'
            specificity_CI_table_from_CI.at[weak_threshold, strong_threshold] = '[' + str(specificity_CI_from_CI_tables['lower_CL'][0]) + ', ' + str(specificity_CI_from_CI_tables['upper_CL'][0]) + ']'
            precision_table_from_CI.at[weak_threshold, strong_threshold]   = stat_scores_from_CI_tables['Precision'][0]
            accuracy_table_from_CI.at[weak_threshold, strong_threshold]    = stat_scores_from_CI_tables['Accuracy'][0]
            MCC_table_from_CI.at[weak_threshold, strong_threshold]         = stat_scores_from_CI_tables['MCC'][0]

            sensitivity_table.at[weak_threshold, strong_threshold] = stat_scores_tables['Sensitivity'][0]
            specificity_table.at[weak_threshold, strong_threshold] = stat_scores_tables['Specificity'][0]
            sensitivity_CI_table.at[weak_threshold, strong_threshold] = '[' + str(sensitivity_CI_tables['lower_CL'][0]) + ', ' + str(sensitivity_CI_tables['upper_CL'][0]) + ']'
            specificity_CI_table.at[weak_threshold, strong_threshold] = '[' + str(specificity_CI_tables['lower_CL'][0]) + ', ' + str(specificity_CI_tables['upper_CL'][0]) + ']'
            precision_table.at[weak_threshold, strong_threshold]   = stat_scores_tables['Precision'][0]
            accuracy_table.at[weak_threshold, strong_threshold]    = stat_scores_tables['Accuracy'][0]
            MCC_table.at[weak_threshold, strong_threshold]         = stat_scores_tables['MCC'][0]

            for signature in all_signatures:
                sensitivity_tables_per_sig_from_CI[signature].at[weak_threshold, strong_threshold] = stat_scores_from_CI_per_sig[signature]['Sensitivity'][0]
                specificity_tables_per_sig_from_CI[signature].at[weak_threshold, strong_threshold] = stat_scores_from_CI_per_sig[signature]['Specificity'][0]
                sensitivity_CI_tables_per_sig_from_CI[signature].at[weak_threshold, strong_threshold] = '[' + str(sensitivity_CI_from_CI_per_sig[signature]['lower_CL'][0]) + ', ' + str(sensitivity_CI_from_CI_per_sig[signature]['upper_CL'][0]) + ']'
                specificity_CI_tables_per_sig_from_CI[signature].at[weak_threshold, strong_threshold] = '[' + str(specificity_CI_from_CI_per_sig[signature]['lower_CL'][0]) + ', ' + str(specificity_CI_from_CI_per_sig[signature]['upper_CL'][0]) + ']'
                precision_tables_per_sig_from_CI[signature].at[weak_threshold, strong_threshold] = stat_scores_from_CI_per_sig[signature]['Precision'][0]
                accuracy_tables_per_sig_from_CI[signature].at[weak_threshold, strong_threshold] = stat_scores_from_CI_per_sig[signature]['Accuracy'][0]
                MCC_tables_per_sig_from_CI[signature].at[weak_threshold, strong_threshold] = stat_scores_from_CI_per_sig[signature]['MCC'][0]
                sensitivity_tables_per_sig[signature].at[weak_threshold, strong_threshold] = stat_scores_per_sig[signature]['Sensitivity'][0]
                specificity_tables_per_sig[signature].at[weak_threshold, strong_threshold] = stat_scores_per_sig[signature]['Specificity'][0]
                sensitivity_CI_tables_per_sig[signature].at[weak_threshold, strong_threshold] = '[' + str(sensitivity_CI_per_sig[signature]['lower_CL'][0]) + ', ' + str(sensitivity_CI_per_sig[signature]['upper_CL'][0]) + ']'
                specificity_CI_tables_per_sig[signature].at[weak_threshold, strong_threshold] = '[' + str(specificity_CI_per_sig[signature]['lower_CL'][0]) + ', ' + str(specificity_CI_per_sig[signature]['upper_CL'][0]) + ']'
                precision_tables_per_sig[signature].at[weak_threshold, strong_threshold] = stat_scores_per_sig[signature]['Precision'][0]
                accuracy_tables_per_sig[signature].at[weak_threshold, strong_threshold] = stat_scores_per_sig[signature]['Accuracy'][0]
                MCC_tables_per_sig[signature].at[weak_threshold, strong_threshold] = stat_scores_per_sig[signature]['MCC'][0]


    # similarity metrics with simple approach
    similarity_tables, similarity_uncertainty_tables, rss_table, chi2_table = measure_similarity_metrics(method)

    # calculate optimal penalties and write to files
    optimal_penalties = pd.DataFrame(dtype=np.float64, index=all_signatures, columns = ['optimal_weak_penalty', 'optimal_strong_penalty'])
    weak_penalty_file = open(output_path + "/optimal_weak_penalty", "w+")
    strong_penalty_file = open(output_path + "/optimal_strong_penalty", "w+")
    # use tables according to CI/no-CI strategy
    if no_CI:
        sensitivity_table_to_use = sensitivity_table
        specificity_table_to_use = specificity_table
        sensitivity_tables_per_sig_to_use = sensitivity_tables_per_sig
        specificity_tables_per_sig_to_use = specificity_tables_per_sig
    else:
        sensitivity_table_to_use = sensitivity_table_from_CI
        specificity_table_to_use = specificity_table_from_CI
        sensitivity_tables_per_sig_to_use = sensitivity_tables_per_sig_from_CI
        specificity_tables_per_sig_to_use = specificity_tables_per_sig_from_CI
    # penalties calculation and filling dataframe
    for signature in all_signatures:
        optimal_weak_penalty, optimal_strong_penalty = calculate_optimal_penalty(sensitivity_tables_per_sig_to_use[signature], specificity_tables_per_sig_to_use[signature], label = signature, metric_to_prioritise=metric_to_prioritise, threshold=metric_threshold)
        optimal_penalties.at[signature, 'optimal_weak_penalty'] = optimal_weak_penalty
        optimal_penalties.at[signature, 'optimal_strong_penalty'] = optimal_strong_penalty
    average_optimal_weak_penalty, average_optimal_strong_penalty = calculate_optimal_penalty(sensitivity_table_to_use, specificity_table_to_use, label = 'average', metric_to_prioritise=metric_to_prioritise, threshold=metric_threshold)
    if signatures_to_prioritise:
        # use with care for multiple signatures
        optimal_weak_penalty = optimal_penalties.loc[signatures_to_prioritise,:]['optimal_weak_penalty'].dropna().max()
        optimal_strong_penalty = optimal_penalties.loc[signatures_to_prioritise,:]['optimal_strong_penalty'].dropna().max()
        # on average calculation for prioritised signatures: not implemented
    else:
        optimal_weak_penalty = optimal_penalties['optimal_weak_penalty'].dropna().max()
        optimal_strong_penalty = optimal_penalties['optimal_strong_penalty'].dropna().max()
        if average:
            optimal_weak_penalty = average_optimal_weak_penalty
            optimal_strong_penalty = average_optimal_strong_penalty
    optimal_penalties.at['average', 'optimal_weak_penalty'] = average_optimal_weak_penalty
    optimal_penalties.at['average', 'optimal_strong_penalty'] = average_optimal_strong_penalty
    optimal_penalties.at['selected', 'optimal_weak_penalty'] = optimal_weak_penalty
    optimal_penalties.at['selected', 'optimal_strong_penalty'] = optimal_strong_penalty
    weak_penalty_file.write("%s\n" % optimal_weak_penalty)
    strong_penalty_file.write("%s\n" % optimal_strong_penalty)

    # summarising metrics for selected penalties
    selected_penalty_metrics = pd.DataFrame(dtype=np.float64, index=all_signatures, columns = ['Sensitivity', 'Specificity', 'Precision', 'Accuracy', 'MCC'])
    selected_penalty_metrics_from_CI = pd.DataFrame(dtype=np.float64, index=all_signatures, columns = ['Sensitivity', 'Specificity', 'Precision', 'Accuracy', 'MCC'])
    for signature in all_signatures:
        selected_penalty_metrics.at[signature, 'Sensitivity'] = np.nan if (optimal_weak_penalty!=optimal_weak_penalty or optimal_strong_penalty!=optimal_strong_penalty) else sensitivity_tables_per_sig[signature].at[optimal_weak_penalty, optimal_strong_penalty]
        selected_penalty_metrics.at[signature, 'Specificity'] = np.nan if (optimal_weak_penalty!=optimal_weak_penalty or optimal_strong_penalty!=optimal_strong_penalty) else specificity_tables_per_sig[signature].at[optimal_weak_penalty, optimal_strong_penalty]
        selected_penalty_metrics.at[signature, 'Precision'] = np.nan if (optimal_weak_penalty!=optimal_weak_penalty or optimal_strong_penalty!=optimal_strong_penalty) else precision_tables_per_sig[signature].at[optimal_weak_penalty, optimal_strong_penalty]
        selected_penalty_metrics.at[signature, 'Accuracy'] = np.nan if (optimal_weak_penalty!=optimal_weak_penalty or optimal_strong_penalty!=optimal_strong_penalty) else accuracy_tables_per_sig[signature].at[optimal_weak_penalty, optimal_strong_penalty]
        selected_penalty_metrics.at[signature, 'MCC'] = np.nan if (optimal_weak_penalty!=optimal_weak_penalty or optimal_strong_penalty!=optimal_strong_penalty) else MCC_tables_per_sig[signature].at[optimal_weak_penalty, optimal_strong_penalty]
        selected_penalty_metrics_from_CI.at[signature, 'Sensitivity'] = np.nan if (optimal_weak_penalty!=optimal_weak_penalty or optimal_strong_penalty!=optimal_strong_penalty) else sensitivity_tables_per_sig_from_CI[signature].at[optimal_weak_penalty, optimal_strong_penalty]
        selected_penalty_metrics_from_CI.at[signature, 'Specificity'] = np.nan if (optimal_weak_penalty!=optimal_weak_penalty or optimal_strong_penalty!=optimal_strong_penalty) else specificity_tables_per_sig_from_CI[signature].at[optimal_weak_penalty, optimal_strong_penalty]
        selected_penalty_metrics_from_CI.at[signature, 'Precision'] = np.nan if (optimal_weak_penalty!=optimal_weak_penalty or optimal_strong_penalty!=optimal_strong_penalty) else precision_tables_per_sig_from_CI[signature].at[optimal_weak_penalty, optimal_strong_penalty]
        selected_penalty_metrics_from_CI.at[signature, 'Accuracy'] = np.nan if (optimal_weak_penalty!=optimal_weak_penalty or optimal_strong_penalty!=optimal_strong_penalty) else accuracy_tables_per_sig_from_CI[signature].at[optimal_weak_penalty, optimal_strong_penalty]
        selected_penalty_metrics_from_CI.at[signature, 'MCC'] = np.nan if (optimal_weak_penalty!=optimal_weak_penalty or optimal_strong_penalty!=optimal_strong_penalty) else MCC_tables_per_sig_from_CI[signature].at[optimal_weak_penalty, optimal_strong_penalty]

    # write dictionaries and dataframes to files
    optimal_penalties.to_csv(output_path + '/optimal_penalties.csv', na_rep='NAN')
    selected_penalty_metrics.to_csv(output_path + '/selected_penalty_metrics.csv', na_rep='NAN')
    selected_penalty_metrics_from_CI.to_csv(output_path + '/selected_penalty_metrics_from_CI.csv', na_rep='NAN')
    sensitivity_table.to_csv(output_path + '/sensitivity_table.csv')
    specificity_table.to_csv(output_path + '/specificity_table.csv')
    sensitivity_CI_table.to_csv(output_path + '/sensitivity_CI_table.csv')
    specificity_CI_table.to_csv(output_path + '/specificity_CI_table.csv')
    precision_table.to_csv(output_path + '/precision_table.csv')
    accuracy_table.to_csv(output_path + '/accuracy_table.csv')
    MCC_table.to_csv(output_path + '/MCC_table.csv')
    rss_table.to_csv(output_path + '/rss_table.csv')
    chi2_table.to_csv(output_path + '/chi2_table.csv')
    write_data_to_JSON(similarity_tables, output_path + '/similarity_tables.json')
    write_data_to_JSON(similarity_uncertainty_tables, output_path + '/similarity_uncertainty_tables.json')

    write_data_to_JSON(sensitivity_tables_per_sig, output_path + '/sensitivity_tables_per_sig.json')
    write_data_to_JSON(specificity_tables_per_sig, output_path + '/specificity_tables_per_sig.json')
    write_data_to_JSON(sensitivity_CI_tables_per_sig, output_path + '/sensitivity_CI_tables_per_sig.json')
    write_data_to_JSON(specificity_CI_tables_per_sig, output_path + '/specificity_CI_tables_per_sig.json')
    write_data_to_JSON(precision_tables_per_sig, output_path + '/precision_tables_per_sig.json')
    write_data_to_JSON(accuracy_tables_per_sig, output_path + '/accuracy_tables_per_sig.json')
    write_data_to_JSON(MCC_tables_per_sig, output_path + '/MCC_tables_per_sig.json')

    sensitivity_table_from_CI.to_csv(output_path + '/sensitivity_table_from_CI.csv')
    specificity_table_from_CI.to_csv(output_path + '/specificity_table_from_CI.csv')
    sensitivity_CI_table_from_CI.to_csv(output_path + '/sensitivity_CI_table_from_CI.csv')
    specificity_CI_table_from_CI.to_csv(output_path + '/specificity_CI_table_from_CI.csv')
    precision_table_from_CI.to_csv(output_path + '/precision_table_from_CI.csv')
    accuracy_table_from_CI.to_csv(output_path + '/accuracy_table_from_CI.csv')
    MCC_table_from_CI.to_csv(output_path + '/MCC_table_from_CI.csv')

    write_data_to_JSON(sensitivity_tables_per_sig_from_CI, output_path + '/sensitivity_tables_per_sig_from_CI.json')
    write_data_to_JSON(specificity_tables_per_sig_from_CI, output_path + '/specificity_tables_per_sig_from_CI.json')
    write_data_to_JSON(sensitivity_CI_tables_per_sig_from_CI, output_path + '/sensitivity_CI_tables_per_sig_from_CI.json')
    write_data_to_JSON(specificity_CI_tables_per_sig_from_CI, output_path + '/specificity_CI_tables_per_sig_from_CI.json')
    write_data_to_JSON(precision_tables_per_sig_from_CI, output_path + '/precision_tables_per_sig_from_CI.json')
    write_data_to_JSON(accuracy_tables_per_sig_from_CI, output_path + '/accuracy_tables_per_sig_from_CI.json')
    write_data_to_JSON(MCC_tables_per_sig_from_CI, output_path + '/MCC_tables_per_sig_from_CI.json')
