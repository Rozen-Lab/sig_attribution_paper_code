import argparse
import os
import copy
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from common_methods import make_folder_if_not_exists, read_data_from_JSON
from ast import literal_eval
matplotlib.use("agg")

def make_lineplot(input_dataframe, xlabel, ylabel, title='', savepath="./lineplot.pdf"):
    global confidence_level
    make_folder_if_not_exists(savepath.rsplit('/', 1)[0])

    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)
    ax.set_title(title, fontsize=14, pad=10)
    ax.set_xlabel(xlabel, size=12)
    ax.set_ylabel(ylabel, size=12)
    ax.tick_params(axis='x', rotation=90)

    for column in input_dataframe.columns:
        if '_CI' in column:
            continue
        elif 'Sensitivity' in column and 'Sensitivity_CI' in input_dataframe.columns:
            errorbars = np.transpose(np.array([literal_eval(CI) for CI in input_dataframe['Sensitivity_CI'].str.replace('nan', '0').to_numpy()]))
            errorbars = np.abs(errorbars-np.array(input_dataframe[column])[None, :])
            plt.errorbar(input_dataframe.index, input_dataframe[column], errorbars, fmt='--o', capsize=1, zorder=10, label=column)
        elif 'Specificity' in column and 'Specificity_CI' in input_dataframe.columns:
            errorbars = np.transpose(np.array([literal_eval(CI) for CI in input_dataframe['Specificity_CI'].str.replace('nan', '0').to_numpy()]))
            errorbars = np.abs(errorbars-np.array(input_dataframe[column])[None, :])
            plt.errorbar(input_dataframe.index, input_dataframe[column], errorbars, fmt='--o', capsize=1, zorder=11, label=column)
        else:
            plt.errorbar(input_dataframe.index, input_dataframe[column], 0, fmt='--o', label=column)

    # Add a grid to the plot
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax.xaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    # Add a red line at the confidence level
    plt.axhline(y=confidence_level, color='r', linestyle='-')

    # legend
    plt.legend(title='Metrics')

    plt.tight_layout()
    plt.savefig(savepath, transparent=True)
    plt.close()

def make_heatmap_plot(input_table, title='', x_label='', y_label='', colormap_label='', savepath="./heatmap_test.pdf", cmap='YlGnBu'):
    table = copy.deepcopy(input_table)
    output_folder = savepath.rsplit('/', 1)[0]
    make_folder_if_not_exists(output_folder)

    if abs(table.values.max()) > 10 or abs(table.values.max()) < 0.01:
        annotation_format = ".1g"
        annotation_text_size = 8
    else:
        annotation_format = ".2g"
        annotation_text_size = 6

    if 'diff' in colormap_label:
        min_matrix = table.values.min()
        max_matrix = table.values.max()
        limit = max(abs(min_matrix), abs(max_matrix))
        cmap = 'bwr'
        annotation_text_size = 8
        ax = sns.heatmap(table, annot=True, fmt=annotation_format, cmap=cmap, vmin=-limit, vmax=limit, linewidths=0.5, xticklabels=True,
                     yticklabels=True, annot_kws={"size": annotation_text_size}, cbar_kws={'label': colormap_label})
    else:
        ax = sns.heatmap(table, annot=True, fmt=annotation_format, cmap=cmap, linewidths=0.5, xticklabels=True,
                     yticklabels=True, annot_kws={"size": annotation_text_size}, cbar_kws={'label': colormap_label})
    ax.set_title(title, fontsize=14, pad=10)

    plt.xlabel(x_label, size=12)
    plt.ylabel(y_label, size=12)

    plt.tight_layout()
    figure = ax.get_figure()
    figure.savefig(savepath)
    plt.close()


def plot_roc_curves(sensitivity, specificity, x_label = '1 - Specificity', y_label = 'Sensitivity', title = 'ROC curves', savepath = 'test_roc.pdf'):
    true_positive_rate = sensitivity
    false_positive_rate = 1 - specificity

    # print(true_positive_rate)
    # print(false_positive_rate)

    # ax = plt.axes()
    curves = plt.plot(false_positive_rate, true_positive_rate)
    dots = plt.scatter(false_positive_rate, true_positive_rate)
    plt.legend(curves, false_positive_rate.index.to_list(), title='Strong threshold')

    plt.plot([0,1], [0,1], ls='--')

    plt.title(title)
    plt.xlabel(x_label, size=12)
    plt.ylabel(y_label, size=12)

    plt.xlim(0,1)
    plt.ylim(0,1)

    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--mutation_type", dest="mutation_type", default='',
                      help="set mutation type (SBS, DBS, ID, SV, CNV)")
    parser.add_argument("-d", "--dataset", dest="dataset_name", default='SIM',
                      help="set the dataset name ('SIM' by default)")
    parser.add_argument("-c", "--context", dest="context", default=96, type=int,
                      help="set SBS context (96, 192, 288, 1536)")
    parser.add_argument("-l", "--confidence_level", dest="confidence_level", default=0.95, type=float,
                      help="specify the confidence level for CL line plotting (default: 0.95)")
    parser.add_argument("-m", "--method", dest="method", default='NNLS',
                      help="set the method name (e.g. 'NNLS')")
    parser.add_argument("-i", "--input_path", dest="input_path", default='output_optimisation/',
                      help="set path to optimisation tables")
    parser.add_argument("-o", "--output_path", dest="output_path", default='efficiency_plots/',
                      help="set path to save output plots")
    parser.add_argument("-s", "--suffix", dest="suffix", default='',
                      help="set a suffix to add to plot names")
    parser.add_argument("-n", "--number", dest="number_of_samples", default=-1, type=int,
                      help="set the number of samples to consider (all by default)")
    parser.add_argument("-W", "--weak_thresholds", nargs='+', dest="weak_thresholds",
                      help="set the list of weak thresholds")
    parser.add_argument("-S", "--strong_thresholds", nargs='+', dest="strong_thresholds",
                      help="set the list of strong thresholds")
    parser.add_argument("--signature_path", dest="signature_tables_path", default='signature_tables/',
                      help="set path to signature tables")
    parser.add_argument("-p", "--signature_prefix", dest="signatures_prefix", default='sigProfiler',
                      help="set prefix in signature filenames (sigProfiler by default)")
    parser.add_argument("-T", "--signature_attribution_thresholds", nargs='+', dest="signature_attribution_thresholds",
                      help="set the list of signature attribution thresholds for confidence probability plots (truth study)")
    args = parser.parse_args()

    mutation_type = args.mutation_type
    dataset = args.dataset_name
    method = args.method
    context = args.context
    confidence_level = args.confidence_level
    signature_tables_path = args.signature_tables_path
    signatures_prefix = args.signatures_prefix
    signature_attribution_thresholds = [float(i)/100 for i in args.signature_attribution_thresholds]

    weak_thresholds = args.weak_thresholds
    strong_thresholds = args.strong_thresholds

    if not weak_thresholds or not strong_thresholds:
        raise ValueError("Please specify the lists of thresholds.")

    if not mutation_type:
        parser.error("Please specify the mutation type using -t option, e.g. add '-t SBS' to the command (DBS, ID).")
    elif mutation_type not in ['SBS','DBS','ID','SV','CNV']:
        raise ValueError("Unknown mutation type: %s. Known types: SBS, DBS, ID, SV, CNV" % mutation_type)

    metrics = ['Cosine', 'Correlation', 'L1', 'L2', 'Chebyshev', 'Jensen-Shannon']

    input_path = args.input_path + '/' + dataset + '/' + mutation_type
    output_path = args.output_path + '/' + dataset + '/' + mutation_type
    number_of_samples = args.number_of_samples

    make_folder_if_not_exists(output_path)

    if mutation_type=='SBS':
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

    # truth_table = pd.read_csv(truth_table_filename, index_col=0)
    all_signatures = list(signatures.columns)
    print('Signatures plotted:', all_signatures)

    # read data from files
    sensitivity_table = pd.read_csv(input_path + '/sensitivity_table.csv', index_col=0)
    specificity_table = pd.read_csv(input_path + '/specificity_table.csv', index_col=0)
    sensitivity_CI_table = pd.read_csv(input_path + '/sensitivity_CI_table.csv', index_col=0)
    specificity_CI_table = pd.read_csv(input_path + '/specificity_CI_table.csv', index_col=0)
    precision_table = pd.read_csv(input_path + '/precision_table.csv', index_col=0)
    accuracy_table = pd.read_csv(input_path + '/accuracy_table.csv', index_col=0)
    MCC_table = pd.read_csv(input_path + '/MCC_table.csv', index_col=0)
    rss_table = pd.read_csv(input_path + '/rss_table.csv', index_col=0)
    chi2_table = pd.read_csv(input_path + '/chi2_table.csv', index_col=0)
    similarity_tables = read_data_from_JSON(input_path + '/similarity_tables.json')
    similarity_uncertainty_tables = read_data_from_JSON(input_path + '/similarity_uncertainty_tables.json')

    sensitivity_tables_per_sig = read_data_from_JSON(input_path + '/sensitivity_tables_per_sig.json')
    specificity_tables_per_sig = read_data_from_JSON(input_path + '/specificity_tables_per_sig.json')
    sensitivity_CI_tables_per_sig = read_data_from_JSON(input_path + '/sensitivity_CI_tables_per_sig.json')
    specificity_CI_tables_per_sig = read_data_from_JSON(input_path + '/specificity_CI_tables_per_sig.json')
    precision_tables_per_sig = read_data_from_JSON(input_path + '/precision_tables_per_sig.json')
    accuracy_tables_per_sig = read_data_from_JSON(input_path + '/accuracy_tables_per_sig.json')
    MCC_tables_per_sig = read_data_from_JSON(input_path + '/MCC_tables_per_sig.json')

    sensitivity_table_from_CI = pd.read_csv(input_path + '/sensitivity_table_from_CI.csv', index_col=0)
    specificity_table_from_CI = pd.read_csv(input_path + '/specificity_table_from_CI.csv', index_col=0)
    sensitivity_CI_table_from_CI = pd.read_csv(input_path + '/sensitivity_CI_table_from_CI.csv', index_col=0)
    specificity_CI_table_from_CI = pd.read_csv(input_path + '/specificity_CI_table_from_CI.csv', index_col=0)
    precision_table_from_CI = pd.read_csv(input_path + '/precision_table_from_CI.csv', index_col=0)
    accuracy_table_from_CI = pd.read_csv(input_path + '/accuracy_table_from_CI.csv', index_col=0)
    MCC_table_from_CI = pd.read_csv(input_path + '/MCC_table_from_CI.csv', index_col=0)

    sensitivity_tables_per_sig_from_CI = read_data_from_JSON(input_path + '/sensitivity_tables_per_sig_from_CI.json')
    specificity_tables_per_sig_from_CI = read_data_from_JSON(input_path + '/specificity_tables_per_sig_from_CI.json')
    sensitivity_CI_tables_per_sig_from_CI = read_data_from_JSON(input_path + '/sensitivity_CI_tables_per_sig_from_CI.json')
    specificity_CI_tables_per_sig_from_CI = read_data_from_JSON(input_path + '/specificity_CI_tables_per_sig_from_CI.json')
    precision_tables_per_sig_from_CI = read_data_from_JSON(input_path + '/precision_tables_per_sig_from_CI.json')
    accuracy_tables_per_sig_from_CI = read_data_from_JSON(input_path + '/accuracy_tables_per_sig_from_CI.json')
    MCC_tables_per_sig_from_CI = read_data_from_JSON(input_path + '/MCC_tables_per_sig_from_CI.json')

    # plot 'ROC' curves using sensitivity/specificity info
    plot_roc_curves(sensitivity_table, specificity_table, savepath=output_path + '/roc_curves_' + str(context) + '_' + method + '.pdf')
    make_folder_if_not_exists(output_path + '/CP_curves')

    signatures_CPs_wrt_thresholds_dict = {}
    for weak_threshold in weak_thresholds:
        for strong_threshold in strong_thresholds:
            input_attributions_folder = args.input_path + '/' + dataset + '_' + str(context) + '_' + method + '_' + weak_threshold + '_' + strong_threshold
            if len(strong_thresholds)>1:
                thresholds_string = weak_threshold + '_' + strong_threshold
            else:
                thresholds_string = weak_threshold
            signatures_CPs = pd.read_csv(input_attributions_folder + '/truth_studies/signatures_CPs_' + mutation_type + '.csv', index_col=0)
            signatures_CPs_dict = read_data_from_JSON(input_attributions_folder + '/truth_studies/signatures_CPs_dict_' + mutation_type + '.json')
            signatures_CPs_wrt_thresholds_dict[thresholds_string] = signatures_CPs_dict

    # condidence probability lines
    for signature in all_signatures:
        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_subplot(111)
        ax.set_title('Signature %s' % signature, fontsize=14, pad=10)
        ax.set_xlabel('Sensitivity threshold', size=12)
        ax.set_ylabel('Confidence probability', size=12)
        ax.tick_params(axis='x', rotation=90)
        ax.grid(True)
        # print(signatures_CPs_wrt_thresholds_dict)

        for key, signature_CP_per_parameters in signatures_CPs_wrt_thresholds_dict.items():
            CPs = [signature_CP_per_parameters[str(threshold)][signature][0] for threshold in signature_attribution_thresholds]
            color = next(ax._get_lines.prop_cycler)['color']
            plt.axhline(y=confidence_level, color='r', linestyle='--')
            plt.scatter(signature_attribution_thresholds, CPs, label=key, color=color, s=10)
            plt.plot(signature_attribution_thresholds, CPs, color=color, linestyle='--')
        plt.legend(loc="best")
        plt.tight_layout()
        plt.savefig(output_path + '/CP_curves/signature_%s_curve.pdf' % signature, transparent=True)
        plt.close()


        # # plot sensitivity threshold vs CP plot for each sig
        # for signature in signatures_to_consider:
        #     signature_attribution_thresholds = sorted(signatures_CPs_dict.keys())
        #     CPs = [signatures_CPs_dict[threshold][signature][0] for threshold in signature_attribution_thresholds]
        #     make_lineplot(signature_attribution_thresholds, CPs, 'Sensitivity threshold', 'Confidence probability', 'Signature %s' % signature, output_folder + '/truth_studies/CP_curves/signature_%s_curve.pdf' % signature)

        # make_boxplot(signatures_CPs, 'Confidence probability (%s)' % mutation_type, 'Signatures', 'CP (%)', savepath=output_folder + '/truth_studies/signature_confidence_probabilities.pdf')
        # make_boxplot(sensitivity_thresholds, 'Sensitivity thresholds (%s)' % mutation_type, 'Signatures', 'Threshold (%)', ylim_zero_to_one = False, savepath=output_folder + '/truth_studies/sensitivity_thresholds.pdf')
        # # make_boxplot(stat_scores_tables, 'Metrics (%s)' % mutation_type, 'Metrics', 'Values', savepath=output_folder + '/sensitivity_specificity_metrics.pdf')
        # make_boxplot(stat_scores_from_CI_tables, 'Metrics (%s)' % mutation_type, 'Metrics', 'Values', savepath=output_folder + '/truth_studies/sensitivity_specificity_from_CI_metrics.pdf')
        # for signature in signatures_to_consider:
        #     # make_boxplot(stat_scores_per_sig[signature], 'Signature %s metrics (%s)' % (signature, mutation_type), 'Metrics', 'Values', savepath=output_folder + '/per_signature_metrics/%s_metrics.pdf' % signature)
        #     make_boxplot(stat_scores_from_CI_per_sig[signature], 'Signature %s metrics (%s)' % (signature, mutation_type), 'Metrics', 'Values', savepath=output_folder + '/truth_studies/per_signature_from_CI_metrics/%s_metrics.pdf' % signature)
        # for score in scores:
        #     make_boxplot(signatures_scores[score], '%s (%s)' % (score, mutation_type), 'Signatures', 'Values', savepath=output_folder + '/truth_studies/signature_%s.pdf' % score)
        # make_boxplot(truth_and_measured_difference, 'Truth and measured difference (%s)' % mutation_type, '', 'Difference', ylim_zero_to_one=False, add_jitter = True, savepath=output_folder + '/truth_studies/truth_and_measured_difference.pdf')

    # assume a single strong threshold to plot more readable line plots
    if len(strong_thresholds)==1:
        # from CIs
        all_metric_dataframes = [sensitivity_table_from_CI.rename(columns={str(strong_thresholds[0]): "Sensitivity"}), sensitivity_CI_table_from_CI.rename(columns={str(strong_thresholds[0]): "Sensitivity_CI"}), specificity_table_from_CI.rename(columns={str(strong_thresholds[0]): "Specificity"}), specificity_CI_table_from_CI.rename(columns={str(strong_thresholds[0]): "Specificity_CI"}), precision_table_from_CI.rename(columns={str(strong_thresholds[0]): "Precision"}), accuracy_table_from_CI.rename(columns={str(strong_thresholds[0]): "Accuracy"}), MCC_table_from_CI.rename(columns={str(strong_thresholds[0]): "MCC"})]
        combined_metrics_dataframe = pd.DataFrame().join(all_metric_dataframes, how="outer")
        make_lineplot(combined_metrics_dataframe, 'Weak penalty', 'Value', title='Metrics from CIs (overall %s)' % (mutation_type), savepath=output_path + "/combined_metrics_from_CIs.pdf")
        # regular
        all_metric_dataframes = [sensitivity_table.rename(columns={str(strong_thresholds[0]): "Sensitivity"}), sensitivity_CI_table.rename(columns={str(strong_thresholds[0]): "Sensitivity_CI"}), specificity_table.rename(columns={str(strong_thresholds[0]): "Specificity"}), specificity_CI_table.rename(columns={str(strong_thresholds[0]): "Specificity_CI"}), precision_table.rename(columns={str(strong_thresholds[0]): "Precision"}), accuracy_table.rename(columns={str(strong_thresholds[0]): "Accuracy"}), MCC_table.rename(columns={str(strong_thresholds[0]): "MCC"})]
        combined_metrics_dataframe = pd.DataFrame().join(all_metric_dataframes, how="outer")
        make_lineplot(combined_metrics_dataframe, 'Weak penalty', 'Value', title='Metrics (overall %s)' % mutation_type, savepath=output_path + "/combined_metrics.pdf")
        for signature in all_signatures:
            # from CIs
            all_metric_dataframes = [sensitivity_tables_per_sig_from_CI[signature].rename(columns={str(strong_thresholds[0]): "Sensitivity"}), sensitivity_CI_tables_per_sig_from_CI[signature].rename(columns={str(strong_thresholds[0]): "Sensitivity_CI"}), specificity_tables_per_sig_from_CI[signature].rename(columns={str(strong_thresholds[0]): "Specificity"}), specificity_CI_tables_per_sig_from_CI[signature].rename(columns={str(strong_thresholds[0]): "Specificity_CI"}), precision_tables_per_sig_from_CI[signature].rename(columns={str(strong_thresholds[0]): "Precision"}), accuracy_tables_per_sig_from_CI[signature].rename(columns={str(strong_thresholds[0]): "Accuracy"}), MCC_tables_per_sig_from_CI[signature].rename(columns={str(strong_thresholds[0]): "MCC"})]
            combined_metrics_dataframe = pd.DataFrame().join(all_metric_dataframes, how="outer")
            make_lineplot(combined_metrics_dataframe, 'Weak penalty', 'Value', title=signature + ' metrics from CIs', savepath=output_path + "/per_signature_from_CI/" + signature + "_metrics_from_CIs.pdf")
            # regular
            all_metric_dataframes = [sensitivity_tables_per_sig[signature].rename(columns={str(strong_thresholds[0]): "Sensitivity"}), sensitivity_CI_tables_per_sig[signature].rename(columns={str(strong_thresholds[0]): "Sensitivity_CI"}), specificity_tables_per_sig[signature].rename(columns={str(strong_thresholds[0]): "Specificity"}), specificity_CI_tables_per_sig[signature].rename(columns={str(strong_thresholds[0]): "Specificity_CI"}), precision_tables_per_sig[signature].rename(columns={str(strong_thresholds[0]): "Precision"}), accuracy_tables_per_sig[signature].rename(columns={str(strong_thresholds[0]): "Accuracy"}), MCC_tables_per_sig[signature].rename(columns={str(strong_thresholds[0]): "MCC"})]
            combined_metrics_dataframe = pd.DataFrame().join(all_metric_dataframes, how="outer")
            make_lineplot(combined_metrics_dataframe, 'Weak penalty', 'Value', title=signature + ' metrics', savepath=output_path + "/per_signature/" + signature + "_metrics.pdf")
    else:
        # without bootstrap
        make_heatmap_plot(sensitivity_table, title='Sensitivity (%i context, %s)' % (context, method.replace('_', ' ')), x_label='Strong threshold',
                          y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/heatmap_' + str(context) + '_' + method + '_sensitivity.pdf')
        make_heatmap_plot(specificity_table, title='Specificity (%i context, %s)' % (context, method.replace('_', ' ')), x_label='Strong threshold',
                          y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/heatmap_' + str(context) + '_' + method + '_specificity.pdf')
        make_heatmap_plot(precision_table, title='Precision (%i context, %s)' % (context, method.replace('_', ' ')), x_label='Strong threshold',
                          y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/heatmap_' + str(context) + '_' + method + '_precision.pdf')
        make_heatmap_plot(accuracy_table, title='Accuracy (%i context, %s)' % (context, method.replace('_', ' ')), x_label='Strong threshold',
                          y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/heatmap_' + str(context) + '_' + method + '_accuracy.pdf')
        make_heatmap_plot(MCC_table, title='Matthews correlation coefficient (%i context, %s)' % (context, method.replace('_', ' ')), x_label='Strong threshold',
                          y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/heatmap_' + str(context) + '_' + method + '_MCC.pdf')

        # using bootstrap, calculations from CIs
        make_heatmap_plot(sensitivity_table_from_CI, title='Sensitivity (%i context, %s)' % (context, method.replace('_', ' ')), x_label='Strong threshold',
                          y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/overall_from_CI/heatmap_' + str(context) + '_' + method + '_sensitivity.pdf')
        make_heatmap_plot(specificity_table_from_CI, title='Specificity (%i context, %s)' % (context, method.replace('_', ' ')), x_label='Strong threshold',
                          y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/overall_from_CI/heatmap_' + str(context) + '_' + method + '_specificity.pdf')
        make_heatmap_plot(precision_table_from_CI, title='Precision (%i context, %s)' % (context, method.replace('_', ' ')), x_label='Strong threshold',
                          y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/overall_from_CI/heatmap_' + str(context) + '_' + method + '_precision.pdf')
        make_heatmap_plot(accuracy_table_from_CI, title='Accuracy (%i context, %s)' % (context, method.replace('_', ' ')), x_label='Strong threshold',
                          y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/overall_from_CI/heatmap_' + str(context) + '_' + method + '_accuracy.pdf')
        make_heatmap_plot(MCC_table_from_CI, title='Matthews correlation coefficient (%i context, %s)' % (context, method.replace('_', ' ')), x_label='Strong threshold',
                          y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/overall_from_CI/heatmap_' + str(context) + '_' + method + '_MCC.pdf')

        for signature in all_signatures:
            make_heatmap_plot(sensitivity_tables_per_sig[signature], title='%s sensitivity (%i context, %s)' % (signature, context, method.replace('_', ' ')), x_label='Strong threshold',
                              y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/per_signature/heatmap_' + signature + '_' + str(context) + '_' + method + '_sensitivity.pdf')
            make_heatmap_plot(specificity_tables_per_sig[signature], title='%s specificity (%i context, %s)' % (signature, context, method.replace('_', ' ')), x_label='Strong threshold',
                              y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/per_signature/heatmap_' + signature + '_' + str(context) + '_' + method + '_specificity.pdf')
            make_heatmap_plot(precision_tables_per_sig[signature], title='%s precision (%i context, %s)' % (signature, context, method.replace('_', ' ')), x_label='Strong threshold',
                              y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/per_signature/heatmap_' + signature + '_' + str(context) + '_' + method + '_precision.pdf')
            make_heatmap_plot(accuracy_tables_per_sig[signature], title='%s accuracy (%i context, %s)' % (signature, context, method.replace('_', ' ')), x_label='Strong threshold',
                              y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/per_signature/heatmap_' + signature + '_' + str(context) + '_' + method + '_accuracy.pdf')
            make_heatmap_plot(MCC_tables_per_sig[signature], title='%s Matthews correlation coefficient (%i context, %s)' % (signature, context, method.replace('_', ' ')), x_label='Strong threshold',
                              y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/per_signature/heatmap_' + signature + '_' + str(context) + '_' + method + '_MCC.pdf')
            # from CI
            make_heatmap_plot(sensitivity_tables_per_sig_from_CI[signature], title='%s sensitivity (%i context, %s)' % (signature, context, method.replace('_', ' ')), x_label='Strong threshold',
                              y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/per_signature_from_CI/heatmap_' + signature + '_' + str(context) + '_' + method + '_sensitivity.pdf')
            make_heatmap_plot(specificity_tables_per_sig_from_CI[signature], title='%s specificity (%i context, %s)' % (signature, context, method.replace('_', ' ')), x_label='Strong threshold',
                              y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/per_signature_from_CI/heatmap_' + signature + '_' + str(context) + '_' + method + '_specificity.pdf')
            make_heatmap_plot(precision_tables_per_sig_from_CI[signature], title='%s precision (%i context, %s)' % (signature, context, method.replace('_', ' ')), x_label='Strong threshold',
                              y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/per_signature_from_CI/heatmap_' + signature + '_' + str(context) + '_' + method + '_precision.pdf')
            make_heatmap_plot(accuracy_tables_per_sig_from_CI[signature], title='%s accuracy (%i context, %s)' % (signature, context, method.replace('_', ' ')), x_label='Strong threshold',
                              y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/per_signature_from_CI/heatmap_' + signature + '_' + str(context) + '_' + method + '_accuracy.pdf')
            make_heatmap_plot(MCC_tables_per_sig_from_CI[signature], title='%s Matthews correlation coefficient (%i context, %s)' % (signature, context, method.replace('_', ' ')), x_label='Strong threshold',
                              y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/per_signature_from_CI/heatmap_' + signature + '_' + str(context) + '_' + method + '_MCC.pdf')


    if 'NNLS' in method:
        make_heatmap_plot(rss_table, title='Mean RSS (%i context, %s)' % (context, method), x_label='Strong threshold',
                          y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/heatmap_' + str(context) + '_' + method + '_rss.pdf')
        make_heatmap_plot(chi2_table, title='Mean Chi2 (%i context, %s)' % (context, method), x_label='Strong threshold',
                          y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/heatmap_' + str(context) + '_' + method + '_chi2.pdf')

    for metric in metrics:
        heatmap_savepath = output_path + '/heatmap_' + str(context) + '_' + method + '_' + metric + '.pdf'
        colormap_label = 'Similarity to truth'
        if args.suffix:
            heatmap_savepath = heatmap_savepath.replace('.pdf', '_%s.pdf' % args.suffix)
        # print(metric, similarity_tables[metric])
        make_heatmap_plot(similarity_tables[metric], title=metric + ' metric (%i context, %s)' % (context, method.replace('_', ' ')), x_label='Strong threshold',
                      y_label='Weak threshold', colormap_label=colormap_label, savepath=heatmap_savepath)

        make_heatmap_plot(similarity_uncertainty_tables[metric], title=metric + ' metric uncertainties (%i context, %s)' % (context, method), x_label='Strong threshold',
                          y_label='Weak threshold', colormap_label='Value', savepath=heatmap_savepath.replace('/heatmap', '/errors/heatmap'))
