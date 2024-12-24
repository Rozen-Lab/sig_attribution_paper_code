from optparse import OptionParser
import os, copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from common_methods import make_folder_if_not_exists, calculate_similarity, calculate_stat_scores

def get_comma_separated_args(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(','))

def make_boxplot(input_table, title, xlabel, ylabel, show_mean = False, savepath = "./boxplot.pdf"):
    output_folder = savepath.rsplit('/',1)[0]
    make_folder_if_not_exists(output_folder)

    table = copy.deepcopy(input_table)
    columns = input_table.columns.to_list()

    fig = plt.figure(figsize=(1.3*len(columns), 4))
    ax = fig.add_subplot(111)
    ax.set_title(title, fontsize=14, pad=10)
    ax.set_xlabel(xlabel, size = 12)
    ax.set_ylabel(ylabel, size = 12)
    ax.set_ylim(0, 1.05)
    plt.axhline(y=1, color='r', linestyle='--')
    table.boxplot(column = columns, ax=ax, grid=False, return_type='axes', sym='.',
                 meanline = show_mean, showmeans = show_mean, boxprops = dict(linewidth=3))
    plt.tight_layout()
    plt.savefig(savepath, transparent=True)
    plt.close()

def make_efficiency_comparison_plot(input_data, title, xlabel, ylabel, savepath = "./efficiencies.pdf"):
    output_folder = savepath.rsplit('/',1)[0]
    make_folder_if_not_exists(output_folder)

    if isinstance(input_data, pd.DataFrame):
        methods = input_data.index.to_list()
        metrics = input_data.columns.to_list()
    else:
        methods = input_data.keys()
        metrics = next(iter(input_data.values())).columns.to_list()

    f, axes = plt.subplots(1, len(metrics), sharey=True, figsize=(len(metrics)*1.3, 4))
    f.suptitle(title, fontsize=12, y=0.999)
    axes[0].set_xlabel(xlabel, size = 12)
    axes[0].set_ylabel(ylabel, size = 12)

    for metric, axis in zip(metrics, axes):
        i = 0
        for method in methods:
            i += 1
            if isinstance(input_data, pd.DataFrame):
                data_to_plot = input_data.loc[method, metric]
            else:
                data_to_plot = input_data[method][metric]
            # print(data_to_plot)
            # print(data_to_plot.mean(), data_to_plot.std())
            if isinstance(data_to_plot, float):
                mean = data_to_plot
                err = 0
            else:
                mean = data_to_plot.mean()
                err = data_to_plot.std()
            if 'unoptimised' in method:
                marker = "v"
                linestyle = '--'
            else:
                marker = "o"
                linestyle = '-'
            errorbar = axis.errorbar(i, mean, yerr=err, fmt=marker, label=method)
            errorbar[-1][0].set_linestyle(linestyle)

        axis.set_ylim(0, 1.05)
        axis.axhline(y=1, color='r', linestyle='--')
        axis.set_xlabel(metric, fontsize=10)
        axis.get_xaxis().set_ticks([])


    handles, labels = axes[0].get_legend_handles_labels()
    f.legend(handles, methods, loc='lower center', framealpha=0.0, borderaxespad=0.1, fancybox=False, shadow=False, ncol=2)

    plt.tight_layout()
    f.subplots_adjust(bottom=0.06 + 0.03*len(methods))
    plt.savefig(savepath, transparent=True)
    plt.close()

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-t", "--mutation_type", dest="mutation_type", default='SBS',
                      help="set mutation type (SBS, DBS, ID, SV, CNV)")
    parser.add_option("-d", "--dataset", dest="dataset_name", default='SIM',
                      help="set the dataset name ('SIM' by default)")
    parser.add_option("-i", "--input_reco_path", dest="input_reco_path", default='output_opt_check/',
                      help="set path to input reconstructed weights tables")
    parser.add_option("-I", "--input_truth_path", dest="input_truth_path", default='input_mutation_tables/',
                      help="set path to input truth (simulated) weights tables")
    parser.add_option("-m", "--methods", type='string', action='callback', callback=get_comma_separated_args,
                      dest = "methods", default = ['NNLS_0.0010_0.0010','NNLS_0.0020_0.0020','NNLS_0.0030_0.0030','NNLS_0.0040_0.0040','NNLS_0.0010_0.0010','NNLS_0.0020_0.0020','NNLS_0.0030_0.0030','NNLS_0.0040_0.0040'],
                      help="set input method names (e.g. -d optimised,unoptimised)")
    parser.add_option("-c", "--contexts", type='string', action='callback', callback=get_comma_separated_args,
                      dest = "contexts", default = ['96','96','96','96','192','192','192','192'],
                      help="set contexts corresponding to the input methods, if running SBS (e.g. -c 96,192)")
    parser.add_option("-o", "--output_path", dest="output_path", default='efficiency_plots/',
                      help="set path to save output plots")
    parser.add_option("-s", "--suffix", dest="suffix", default='',
                      help="set a suffix to add to plot names")
    parser.add_option("-n", "--number", dest="number_of_samples", default=-1, type='int',
                      help="set the number of samples to consider (all by default)")
    (options, args) = parser.parse_args()

    mutation_type = options.mutation_type
    dataset = options.dataset_name
    input_reco_path = options.input_reco_path
    input_truth_path = options.input_truth_path + '/' + dataset
    methods = options.methods
    contexts = options.contexts

    metrics = ['Cosine', 'Correlation', 'L1', 'L2', 'Chebyshev', 'Jensen-Shannon']
    scores = ['Sensitivity', 'Specificity', 'Precision', 'Accuracy', 'F1', 'MCC']

    if len(methods)!=len(contexts) and mutation_type=='SBS':
        parser.error("Please verify that provided methods and contexts (-m, -c options) are compatible.")
    if mutation_type!='SBS' and set(methods)!=set(['optimised','unoptimised']):
        # override context:
        print("Warning: overriding methods with ['unoptimised','optimised'] for %s mutation type" % mutation_type)
        methods = ['unoptimised','optimised']

    if mutation_type not in ['SBS','DBS','ID','SV','CNV']:
        raise ValueError("Unknown mutation type: %s. Known types: SBS, DBS, ID, SV, CNV" % mutation_type)

    output_path = options.output_path
    number_of_samples = options.number_of_samples

    make_folder_if_not_exists(output_path)

    if not input_reco_path:
        parser.error("Please specify the input path for reconstructed weights tables using -i option.")
    if not input_truth_path:
        parser.error("Please specify the input truth (simulated) weights table using -I option.")

    if mutation_type=='SBS':
        reco_table_filenames = [input_reco_path + '/' + dataset + '_' + context + '_' + method + '/output_%s_SBS_weights_table.csv' % dataset
                                for context, method in zip(contexts, methods)]
        truth_table_filenames = [input_truth_path + '/WGS_' + dataset + '.' + context + '.weights.csv' for context in contexts]
    elif mutation_type=='DBS':
        reco_table_filenames = [input_reco_path + '/' + dataset + '_' + method + '/output_%s_DBS_weights_table.csv' % dataset for method in methods]
        truth_table_filenames = [input_truth_path + '/WGS_' + dataset + '.dinucs.weights.csv'] * len(methods)
    elif mutation_type=='ID':
        reco_table_filenames = [input_reco_path + '/' + dataset + '_' + method + '/output_%s_ID_weights_table.csv' % dataset for method in methods]
        truth_table_filenames = [input_truth_path + '/WGS_' + dataset + '.indels.weights.csv'] * len(methods)
    else: # SV and CNV
        reco_table_filenames = [input_reco_path + '/' + dataset + '_' + method + '/output_%s_%s_weights_table.csv' % (dataset, mutation_type) for method in methods]
        truth_table_filenames = [input_truth_path + '/WGS_' + dataset + '.' + mutation_type + '.weights.csv'] * len(methods)

    # add contexts to method names
    if mutation_type=='SBS':
        methods = ['_'.join([context, method]) for context, method in zip(contexts, methods)]

    similarity_tables = {}
    stat_scores_tables = pd.DataFrame(index = methods, columns = scores)
    stat_scores_per_sig = {}
    for method, reco_table_filename, truth_table_filename in zip(methods,reco_table_filenames,truth_table_filenames):
        reco_table = pd.read_csv(reco_table_filename, index_col=0)
        truth_table = pd.read_csv(truth_table_filename, index_col=0)
        similarity_tables[method] = pd.DataFrame()

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
        all_signatures = list(truth_table.columns)

        # initialise stat scores dataframes per signature
        for signature in all_signatures:
            if not signature in stat_scores_per_sig.keys():
                stat_scores_per_sig[signature] = pd.DataFrame(index = methods, columns=scores, dtype=float)

        if set(reco_table.columns.to_list()) != set(truth_table.columns.to_list()):
            # print(reco_table.columns.to_list())
            # print(truth_table.columns.to_list())
            raise ValueError("Incompatible truth and reco tables: check signature columns")
        if set(reco_table.index) != set(truth_table.index):
            # print(reco_table.index.to_list())
            # print(truth_table.index.to_list())
            raise ValueError("Incompatible truth and reco tables: check indexing (samples)")

        if number_of_samples == -1 or number_of_samples>len(reco_table):
            number_of_samples = len(reco_table)

        # calculate similarities
        for metric in metrics:
            for i in range(number_of_samples):
                reco_sample = reco_table.iloc[i].values
                truth_sample = truth_table.iloc[i].values
                if not reco_sample.any() and not truth_sample.any():
                    # similarity (e.g. cosine or JS) is undefined for zero samples, so skipping
                    continue
                similarity_tables[method].loc[i, metric] = calculate_similarity(reco_sample, truth_sample, metric)

        # calculate stat scores
        print('Method:', method)
        stat_scores_tables.loc[method, :] = calculate_stat_scores(all_signatures, reco_table, truth_table, number_of_samples)
        for signature in all_signatures:
            stat_scores_per_sig[signature].loc[method, :] = calculate_stat_scores([signature], reco_table, truth_table, number_of_samples)

        boxplot_savepath = output_path + '/box_plots/' + dataset + '_' + mutation_type + '_' + method + '.pdf'
        if options.suffix:
            boxplot_savepath = boxplot_savepath.replace('box_plots','box_plots_%s' % options.suffix)
        make_boxplot(similarity_tables[method], 'Signature attribution efficiency (%s)' % method,
                        'Metrics', 'Similarities', savepath = boxplot_savepath)

    comparison_savepath = output_path + '/' + dataset + '_' + mutation_type + '_comparison_efficiencies.pdf'
    if options.suffix:
        comparison_savepath = comparison_savepath.replace('.pdf','_%s.pdf' % options.suffix)
    make_efficiency_comparison_plot(similarity_tables, 'Signature attribution efficiencies (%s)' % mutation_type, 'Metrics', 'Similarities', savepath = comparison_savepath)
    make_efficiency_comparison_plot(stat_scores_tables, 'Metrics (%s)' % mutation_type, 'Metrics', 'Values', savepath = comparison_savepath.replace('efficiencies.pdf','overall_metrics.pdf'))
    for signature in all_signatures:
        make_efficiency_comparison_plot(stat_scores_per_sig[signature], 'Signature %s metrics' % (signature), 'Metrics', 'Values', savepath = comparison_savepath.replace('efficiencies.pdf', signature + '_metrics.pdf'))
