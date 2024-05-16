import argparse
import os
import copy
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from common_methods import make_folder_if_not_exists, calculate_similarity, calculate_stat_scores
matplotlib.use("agg")

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

def measure_metrics(method):
    global number_of_samples, input_truth_path, input_reco_path
    global dataset, context, metrics, weak_thresholds, strong_thresholds
    # initialise dicts and arrays
    similarity_tables = {}
    similarity_uncertainty_tables = {}

    for metric in metrics:
        similarity_tables[metric] = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
        similarity_uncertainty_tables[metric] = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)

    sensitivity_table = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
    specificity_table = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
    MCC_table = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)

    rss_table = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)
    chi2_table = pd.DataFrame(dtype=np.float64, columns=strong_thresholds, index=weak_thresholds)

    truth_table_filename = input_truth_path + '/WGS_' + dataset + '.' + str(context) + '.weights.csv'
    truth_table = pd.read_csv(truth_table_filename, index_col=0)
    all_signatures = list(truth_table.columns)

    for weak_threshold in weak_thresholds:
        for strong_threshold in strong_thresholds:
            reco_table_filename = input_reco_path + '/' + dataset + '_' + \
                str(context) + '_' + method + '_' + weak_threshold + \
                '_' + strong_threshold + '/output_%s_SBS_weights_table.csv' % dataset
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

            if number_of_samples == -1 or number_of_samples > len(reco_table):
                number_of_samples = len(reco_table)

            # RSS and Chi2 for NNLS from separate stat info table
            if 'NNLS' in method:
                stat_table_filename = input_reco_path + '/' + dataset + '_' + \
                    str(context) + '_' + method + '_' + weak_threshold + \
                    '_' + strong_threshold + '/output_%s_SBS_stat_info.csv' % dataset
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
                    if np.isnan(similarity):
                        print('Warning: Nan similarity in metric %s, thresholds %s/%s' % (metric, weak_threshold, strong_threshold), 'sample:',i,reco_sample,truth_sample)
                    similarities.append(similarity)

                similarity_tables[metric].loc[weak_threshold, strong_threshold] = np.mean(similarities)
                similarity_uncertainty_tables[metric].loc[weak_threshold, strong_threshold] = np.std(
                    similarities)

            # calculate sensitivities and specificities
            sensitivity, specificity, _, _, _, MCC = calculate_stat_scores(all_signatures, reco_table, truth_table, number_of_samples)
            sensitivity_table.loc[weak_threshold, strong_threshold] = sensitivity
            specificity_table.loc[weak_threshold, strong_threshold] = specificity
            MCC_table.loc[weak_threshold, strong_threshold] = MCC

    return sensitivity_table, specificity_table, MCC_table, similarity_tables, similarity_uncertainty_tables, rss_table, chi2_table


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--mutation_type", dest="mutation_type", default='SBS',
                      help="set mutation type (SBS, DBS, ID)")
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
    args = parser.parse_args()

    mutation_type = args.mutation_type
    dataset = args.dataset_name
    input_reco_path = args.input_reco_path
    input_truth_path = args.input_truth_path
    method = args.method
    context = args.context

    weak_thresholds = args.weak_thresholds
    strong_thresholds = args.strong_thresholds
    
    if not weak_thresholds or not strong_thresholds:
        raise ValueError("Please specify the lists of thresholds.")
    ## ultrafine range (adjust according to input weights tables)
    # weak_thresholds = ['%.4f' % (0.0010 * (i + 1)) for i in range(10)]
    # strong_thresholds = ['%.4f' % (0.0010 * (i + 1)) for i in range(10)]
    # weak_thresholds = ['%.3f' % (0.001 * (i + 1)) for i in range(5)]
    # strong_thresholds = ['%.3f' % (0.001 * (i + 1)) for i in range(5)]

    metrics = ['Cosine', 'Correlation', 'L1', 'L2', 'Chebyshev', 'Jensen-Shannon']

    if mutation_type != 'SBS':
        raise ValueError("Unsupported mutation type: %s. Please use SBS." % mutation_type)

    output_path = args.output_path + '/' + dataset
    number_of_samples = args.number_of_samples

    make_folder_if_not_exists(output_path)

    if not method:
        parser.error("Please specify the method name (e.g. 'NNLS', 'SigProfiler') using -m option.")
    if not input_reco_path:
        parser.error("Please specify the input path for reconstructed weights tables using -i option.")
    if not input_truth_path:
        parser.error("Please specify the input truth (simulated) weights table using -I option.")

    sensitivity_table, specificity_table, MCC_table, similarity_tables, similarity_uncertainty_tables, rss_table, chi2_table = measure_metrics(method)

    plot_roc_curves(sensitivity_table, specificity_table, savepath=output_path + '/roc_curves_' + str(context) + '_' + method + '.pdf')

    if 'NNLS' in method:
        make_heatmap_plot(rss_table, title='Mean RSS (%i context, %s)' % (context, method), x_label='Strong threshold',
                          y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/heatmap_' + str(context) + '_' + method + '_rss.pdf')
        make_heatmap_plot(chi2_table, title='Mean Chi2 (%i context, %s)' % (context, method), x_label='Strong threshold',
                          y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/heatmap_' + str(context) + '_' + method + '_chi2.pdf')

    make_heatmap_plot(sensitivity_table, title='Sensitivity (%i context, %s)' % (context, method.replace('_', ' ')), x_label='Strong threshold',
                      y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/heatmap_' + str(context) + '_' + method + '_sensitivity.pdf')
    make_heatmap_plot(specificity_table, title='Specificity (%i context, %s)' % (context, method.replace('_', ' ')), x_label='Strong threshold',
                      y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/heatmap_' + str(context) + '_' + method + '_specificity.pdf')
    make_heatmap_plot(MCC_table, title='Matthews correlation coefficient (%i context, %s)' % (context, method.replace('_', ' ')), x_label='Strong threshold',
                      y_label='Weak threshold', colormap_label='Value', savepath=output_path + '/heatmap_' + str(context) + '_' + method + '_MCC.pdf')

    for metric in metrics:
        heatmap_savepath = output_path + '/heatmap_' + str(context) + '_' + method + '_' + metric + '.pdf'
        colormap_label = 'Similarity to truth'
        if args.suffix:
            heatmap_savepath = heatmap_savepath.replace('.pdf', '_%s.pdf' % args.suffix)
        # print(metric, similarity_tables[metric])
        make_heatmap_plot(similarity_tables[metric], title=metric + ' metric (%i context, %s)' % (context, method.replace('_', ' ')), x_label='Strong threshold',
                      y_label='Weak threshold', colormap_label=colormap_label, savepath=heatmap_savepath)
        make_heatmap_plot(similarity_tables[metric], title=metric + ' metric (%i context, %s)' % (context, method.replace('_', ' ')), x_label='Strong threshold',
                      y_label='Weak threshold', colormap_label=colormap_label, savepath=heatmap_savepath)

        make_heatmap_plot(similarity_uncertainty_tables[metric], title=metric + ' metric uncertainties (%i context, %s)' % (context, method), x_label='Strong threshold',
                          y_label='Weak threshold', colormap_label='Value', savepath=heatmap_savepath.replace('/heatmap', '/errors/heatmap'))
