""" plot_bootstrap_attributions.py
Module to plot attributions with CIs on per-sample and per-signature basis.
"""

from argparse import ArgumentParser
import copy
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from common_methods import make_folder_if_not_exists, read_data_from_JSON
plt.rc('text', usetex=False)

def plot_bootstrap_attribution_histograms(input_table, centrals, title, savepath='./hist.pdf'):
    table = copy.deepcopy(input_table)
    output_folder = savepath.rsplit('/',1)[0]
    make_folder_if_not_exists(output_folder)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    if max(centrals.values())<=1:
        bins = [0.01*k for k in range(101)]
    else:
        bins = None

    for column in table.columns:
        color = next(ax._get_lines.prop_cycler)['color']
        plt.hist(table[column].to_list(), label = column, bins = bins, histtype='step', color=color, density=True)
        plt.axvline(x=centrals[column], color=color, ls = '--')

    legend = plt.legend()
    plt.gca().add_artist(legend)

    ax.set_title(title, fontsize=14, pad=10)
    ax.set_xlabel('Attribution', fontsize=12)
    ax.set_ylabel('Number of samples (normalised)', fontsize=12)
    if max(centrals.values())<=1:
        ax.set_xticks([0.1*k for k in range(11)])

    plt.tight_layout()
    plt.savefig(savepath, transparent=True)
    plt.close()

def make_bootstrap_attribution_boxplot(input_table, centrals = None, truth = None, show_mean = False, prefix = 'Bootstrap', method = 'NNLS', title='', x_label = '', y_label = '', savepath='./boxplot.pdf'):
    table = copy.deepcopy(input_table)
    output_folder = savepath.rsplit('/',1)[0]
    make_folder_if_not_exists(output_folder)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(x_label, size = 12)
    ax.set_ylabel(y_label, size = 12)
    ax.set_title(title, fontsize=14, pad=10)

    sns.set(style="whitegrid")
    ax = sns.boxplot(data=table, palette="colorblind", whis=[2.5, 97.5], showfliers=False, meanline = show_mean, showmeans = show_mean)

    # add counts to labels
    # sig_labels = [t.get_text() for t in ax.get_xticklabels()]
    # sig_labels = [sig + ' (N = %i)' % len(table[sig].dropna()) for sig in sig_labels]
    # ax.set_xticklabels(sig_labels)
    ax.tick_params(axis='x', rotation=90, labelsize=8)
    column_range = range(len(table.columns))

    if centrals or truth:
        # for id, name in zip(column_range, table.columns.to_list()):
        for id in column_range:
            if centrals:
                # plt.plot([column_range[id]], [centrals[name]], marker="o", linestyle="None", markersize=3, color='red', label = method)
                plt.plot([column_range[id]], [list(centrals.values())[id]], marker="o", linestyle="None", markersize=3, color='red', label = method)
            if truth:
                # plt.plot([column_range[id]], [truth[name]], marker="o", linestyle="None", markersize=3, color='green', label = "Truth")
                plt.plot([column_range[id]], [list(truth.values())[id]], marker="o", linestyle="None", markersize=3, color='green', label = "Truth")

    # Add a boxplot patch for legend
    boxplot_patch = mpatches.Patch(color='gray', lw=2, label=prefix + ' ' + method)
    handles, labels = ax.get_legend_handles_labels()
    if truth and centrals:
        handles = handles[:2]
    else:
        handles = handles[:1]

    handles.append(boxplot_patch)
    plt.legend(handles=handles)

    # adjust figure size for long data
    if len(column_range)>100:
        fig = ax.get_figure()
        fig.set_size_inches(20, 5)

    # if relative:
        # ax.set_ylim(0, 1)

    # Add a horizontal grid to the plot
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

    plt.tight_layout()
    plt.savefig(savepath, transparent=True)
    plt.close()

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_attributions_folder", dest="input_attributions_folder", default='output_tables/',
                        help="set path to NNLS output data")
    parser.add_argument("-I", "--input_mutations_folder", dest="input_mutations_folder", default='input_mutation_tables/',
                        help="set path to datasets with input mutation tables")
    parser.add_argument("-S", "--signature_path", dest="signature_tables_path", default='signature_tables/',
                        help="set path to signature tables")
    parser.add_argument("-p", "--signature_prefix", dest="signatures_prefix", default='sigProfiler',
                        help="set prefix in signature filenames (sigProfiler by default)")
    parser.add_argument("-d", "--dataset", dest="dataset_name", default='SIM',
                        help="set the dataset name (e.g. SIM)")
    parser.add_argument("-o", "--output_folder", dest="output_folder", default='plots/',
                        help="set path to save plots")
    parser.add_argument("-t", "--mutation_type", dest="mutation_type", default='',
                        help="set mutation type (SBS, DBS, ID, SV, CNV)")
    parser.add_argument("-c", "--context", dest="context", default=192, type=int,
                        help="set SBS context (96, 192)")
    parser.add_argument("-a", "--plot_absolute_numbers", dest="abs_numbers", action="store_true",
                        help="show absolute numbers of mutations")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                        help="print additional information for debugging")
    parser.add_argument("-n", "--number_of_b_samples", dest="number_of_b_samples", default=1000, type=int,
                        help="Number of bootstrapped samples (1000 by default)")
    parser.add_argument("-N", "--number_of_samples", dest="number_of_samples", default=-1, type=int,
                        help="limit the number of samples to analyse (all by default)")

    args = parser.parse_args()

    dataset_name = args.dataset_name
    mutation_type = args.mutation_type
    context = args.context
    number_of_b_samples = args.number_of_b_samples
    signature_tables_path = args.signature_tables_path
    signatures_prefix = args.signatures_prefix
    input_mutations_folder = args.input_mutations_folder
    input_attributions_folder = args.input_attributions_folder + '/' + dataset_name + '/'
    output_folder = args.output_folder + '/' + dataset_name + '/' + mutation_type + '/bootstrap_plots/'
    make_folder_if_not_exists(output_folder)

    if not mutation_type:
        parser.error("Please specify the mutation type using -t option, e.g. add '-t SBS' to the command (DBS, ID).")
    elif mutation_type not in ['SBS', 'DBS', 'ID', 'SV', 'CNV']:
        raise ValueError("Unknown mutation type: %s. Known types: SBS, DBS, ID, SV, CNV" % mutation_type)

    print("*"*50)
    print("Making plots for %s mutation type, %s dataset" % (mutation_type, dataset_name))

    central_attribution_table_abs = pd.read_csv(input_attributions_folder + '/output_%s_%s_mutations_table.csv' % (dataset_name, mutation_type), index_col=0)
    central_attribution_table_weights = pd.read_csv(input_attributions_folder + '/output_%s_%s_weights_table.csv' % (dataset_name, mutation_type), index_col=0)
    central_stat_table = pd.read_csv(input_attributions_folder + '/output_%s_%s_stat_info.csv' % (dataset_name, mutation_type), index_col=0)

    if mutation_type=='SBS':
        if context==96:
            input_mutations = pd.read_csv('%s/%s/WGS_%s.%i.csv' % (input_mutations_folder, dataset_name, dataset_name, context), index_col=[0,1])
        elif context in [192, 288]:
            input_mutations = pd.read_csv('%s/%s/WGS_%s.%i.csv' % (input_mutations_folder, dataset_name, dataset_name, context), index_col=[0,1,2])
        elif context==1536:
            input_mutations = pd.read_csv('%s/%s/WGS_%s.%i.csv' % (input_mutations_folder, dataset_name, dataset_name, context), index_col=0)
        else:
            raise ValueError("Context %i is not supported." % context)
    elif mutation_type=='DBS':
        input_mutations = pd.read_csv('%s/%s/WGS_%s.dinucs.csv' % (input_mutations_folder, dataset_name, dataset_name), index_col=0)
    elif mutation_type=='ID':
        input_mutations = pd.read_csv('%s/%s/WGS_%s.indels.csv' % (input_mutations_folder, dataset_name, dataset_name), index_col=0)
    else: # SV and CNV
        input_mutations = pd.read_csv('%s/%s/WGS_%s.%s.csv' % (input_mutations_folder, dataset_name, dataset_name, mutation_type), index_col=0)

    if 'SIM' in dataset_name:
        if mutation_type=='SBS':
            truth_attribution_table = pd.read_csv(input_attributions_folder + '/WGS_%s.%i.weights.csv' % (dataset_name, context), index_col=0)
        elif mutation_type=='DBS':
            truth_attribution_table = pd.read_csv(input_attributions_folder + '/WGS_%s.dinucs.weights.csv' % dataset_name, index_col=0)
        elif mutation_type=='ID':
            truth_attribution_table = pd.read_csv(input_attributions_folder + '/WGS_%s.indels.weights.csv' % dataset_name, index_col=0)
        else: # SV and CNV
            truth_attribution_table = pd.read_csv(input_attributions_folder + '/WGS_%s.%s.weights.csv' % (dataset_name, mutation_type), index_col=0)

    if args.abs_numbers:
        if 'SIM' in dataset_name:
            # multiply truth table by mutational burden from input sample:
            truth_attribution_table = truth_attribution_table.mul(np.asarray(input_mutations.sum(axis=0)), axis=0)
        central_attribution_table = central_attribution_table_abs
        colormap_label = 'Absolute mutations number'
        input_filename_suffix = "_%s_%s_bootstrap_output_abs_mutations.json" % (dataset_name, mutation_type)
        savepath_filename = 'bootstrap_plot_abs_mutations'
    else:
        central_attribution_table = central_attribution_table_weights
        colormap_label = 'Relative contribution'
        input_filename_suffix = "_%s_%s_bootstrap_output_weights.json" % (dataset_name, mutation_type)
        savepath_filename = 'bootstrap_plot_weights'

    # limit the number of samples to analyse (if specified by -N option)
    if args.number_of_samples!=-1:
        central_attribution_table = central_attribution_table.head(args.number_of_samples)
    # convert columns and index to str
    central_attribution_table.columns = central_attribution_table.columns.astype(str)
    central_attribution_table.index = central_attribution_table.index.astype(str)

    signatures_to_consider = list(central_attribution_table_abs.columns)

    if 'SIM' in dataset_name:
        truth_attribution_table.columns = truth_attribution_table.columns.astype(str)
        truth_attribution_table.index = truth_attribution_table.index.astype(str)
        for signature in signatures_to_consider:
            # if signature is not modelled, add it as zeros in truth table:
            if signature not in truth_attribution_table.columns:
                truth_attribution_table[signature] = 0

    main_title = dataset_name.replace('_','/') + ' data, ' + mutation_type + ' mutation type'

    # samples and metrics to consider
    samples = central_attribution_table.index.to_list()
    stat_metrics = central_stat_table.columns.to_list()

    # read dictionaries from JSON files produced by make_bootstrap_tables.py script
    stat_metrics_dict = read_data_from_JSON(input_attributions_folder + "stat_metrics" + input_filename_suffix)
    attributions_per_sample_dict = read_data_from_JSON(input_attributions_folder + "attributions_per_sample" + input_filename_suffix)
    attributions_per_signature_dict = read_data_from_JSON(input_attributions_folder + "attributions_per_signature" + input_filename_suffix)
    # # pd.to_json function does not fully support MultiIndex (present in signature tables), hence commented out for now
    # mutation_spectra_dict = read_data_from_JSON(input_attributions_folder + "mutation_spectra" + input_filename_suffix, set_index=[0, 1])

    # make a boxplot with central attributions
    make_bootstrap_attribution_boxplot(central_attribution_table,
        title = main_title + ': NNLS attributions',
        prefix = 'Optimised',
        y_label = 'Signature attribution',
        savepath = output_folder + '/' + savepath_filename.replace('bootstrap_plot','all_attributions') + '.pdf')

    # if applicable, make a boxplot with truth attributions
    if 'SIM' in dataset_name:
        make_bootstrap_attribution_boxplot(truth_attribution_table,
            title = main_title + ': truth attributions',
            prefix = 'Truth',
            method = '',
            y_label = 'Signature attribution',
            savepath = output_folder + '/' + savepath_filename.replace('bootstrap_plot', 'all_true_values') + '.pdf')

    # make bootstrap attribution plots for each sample
    for sample in samples:
        make_bootstrap_attribution_boxplot(attributions_per_sample_dict[sample],
            centrals = central_attribution_table.loc[sample].to_dict(),
            title = main_title + ': sample ' + sample.replace("_", "-"),
            truth = truth_attribution_table.loc[sample].to_dict() if 'SIM' in dataset_name else None,
            y_label = 'Signature attribution',
            savepath = output_folder + '/sample_based/' + savepath_filename + '_' + sample + '.pdf')
        # # pd.to_json function does not fully support MultiIndex (present in signature tables), hence commented out for now
        # make_bootstrap_attribution_boxplot(mutation_spectra_dict[sample],
        #     truth = input_mutations[sample].to_dict(),
        #     title = main_title + ': sample ' + sample.replace("_", "-"),
        #     method = "NNLS (recalculated)",
        #     y_label = 'Mutation count',
        #     savepath = output_folder + '/sample_based/spectra/bootstrap_plot_' + sample + '_spectrum.pdf')
        for signature in signatures_to_consider:
            # do not plot signatures with both central (or truth) and mean bootstrap attributions of zero:
            if 'SIM' in dataset_name:
                if truth_attribution_table.loc[sample, signature] == 0 and np.mean(attributions_per_sample_dict[sample][signature] == 0):
                    attributions_per_sample_dict[sample].drop(signature, axis=1, inplace=True)
            elif central_attribution_table.loc[sample, signature] == 0 and np.mean(attributions_per_sample_dict[sample][signature] == 0):
                attributions_per_sample_dict[sample].drop(signature, axis=1, inplace=True)
        # make sample-based plots with signature attributions as histograms
        plot_bootstrap_attribution_histograms(attributions_per_sample_dict[sample],
            centrals = central_attribution_table.loc[sample].to_dict(),
            title = main_title + ': sample ' + sample.replace("_", "-"),
            savepath = output_folder + '/sample_based/histograms/' + savepath_filename + '_' + sample + '_dist.pdf')

    # make bootstrap attribution plots for each signature
    for signature in signatures_to_consider:
        if 'SIM' in dataset_name:
            if np.mean(truth_attribution_table[signature].values)==0 and np.mean(attributions_per_signature_dict[signature].values)==0:
                # skip zero signatures in both truth and measurement
                continue
        make_bootstrap_attribution_boxplot(attributions_per_signature_dict[signature],
            centrals = central_attribution_table[signature].to_dict(),
            title = main_title + ': signature ' + signature.replace("_", "-"),
            truth = truth_attribution_table[signature].to_dict() if 'SIM' in dataset_name else None,
            y_label = 'Signature attribution',
            savepath = output_folder + '/signature_based/' + savepath_filename + '_' + signature + '.pdf')

    # make bootstrap plots for each statistical metric available
    for metric in stat_metrics:
        make_bootstrap_attribution_boxplot(stat_metrics_dict[metric],
            centrals = central_stat_table[metric].to_dict(),
            title = main_title + ': %s' % (metric),
            x_label = 'Samples',
            y_label = metric,
            method = metric, savepath = output_folder + '/' + metric + '_bootstrap_plot.pdf')
